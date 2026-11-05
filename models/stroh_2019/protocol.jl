# Protocol for stroh_2019: dosing, simulation settings, targets, tolerances, fixes, benchmark setups.
using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq

# ---- targets (frozen 2026-11-04 09:48 local, before any simulation was compared; see notes.md)
# Grade-bearing: Figure 3a-c (monkey calibration lines, 10 curves) and Figure 4a-e (human projections, 25 curves).
# Not graded: observed monkey points, Figure 5.
const TOLERANCE = Tolerance()                       # default rms 0.015, p95 0.03 of the axis box
const TOLERANCES = Dict{String, Tolerance}()        # no per-curve overrides

# ---- fixes (provenance ledger 2026-11-05). Every input the paper does not state, or states only as a bound, is a Fix
# in the run that uses it. The as-published run needs PUBLISHED_FILLS to run at all (without them `simulate_curve`
# throws NotRunnable); the fixed run keeps them and adds `unbound_readout`, the only change between the two runs.
const PRE = " Stage: listed in notes.md under 'Modelling decisions made before comparing'."
const PUBLISHED_FILLS = [
    Fix("molecular_weight_150kDa", :missing,
        "mg/kg doses are converted to nM with MW = 150 kDa (model.jl MW).";
        evidence = "The paper gives no molecular weight. 150 kDa is the standard IgG1 mass (Pb-Tx adds a prodomain, so the true mass is somewhat higher; value is a standard constant, not read from the figure). Consistency only, not the source: 4.5 mg/kg x 70 kg / 150 kDa / 2.6 L = 808 nM." * PRE,
        basis = :standard, stage = :before_comparison),
    Fix("human_V2_2.6L", :missing,
        "Human peripheral volume V2 = 2.6 L (not in Table 2).";
        evidence = "Derived from the paper's own tables: Table 1 monkey receptor concentration 9.0e-5 nmol/s / (1e-4 1/s * 0.1 L) = 9 nM; Results state normal-tissue receptor concentration is assumed equal in human; Table 2 ksynR peripheral 2.3e-3 nmol/s / (1e-4 * 9 nM) = 2.56 L ~ 2.6 L (= V1)." * PRE,
        basis = :paper, stage = :before_comparison, figures = ["fig4"]),
    Fix("kcleave_monkey_3e-7", :ambiguity,
        "Monkey kcleave for substrates S1 and S2 = 3e-7 1/s in plasma and periphery (Table 1 prints '<3e-7'); 0 for M1,0 (no substrate).";
        evidence = "Results: '3e-7 1/second represented an upper bound on this value. This value of 3e-7 1/second was carried forward as the likely maximum value for k cleave' and, for Figure 3, 'with k cleave held fixed to the upper bound of 3e-7 1/second'. The paper names this point." * PRE,
        basis = :interpretation, stage = :before_comparison, figures = ["fig3"]),
    Fix("kcleave_human_3e-7", :ambiguity,
        "Human kcleave in plasma and periphery = 3e-7 1/s (Table 2 prints '<3e-7').";
        evidence = "Table 2: '<3e-7, Assumed to be the same as monkey for central and peripheral compartments'; Results: kcleave 'assumed to have an upper bound of 3e-7'. The monkey value used for Figure 3 is 3e-7 (kcleave_monkey_3e-7). The value used for Figure 4 is not stated; 1e-7 would match it (report.md, Sensitivity), which is not used." * PRE,
        basis = :interpretation, stage = :before_comparison, figures = ["fig4"]),
    Fix("tumour_kcleave_10fold", :ambiguity,
        "Human tumour kcleave = 10 x the plasma/periphery kcleave, in all Figure 4 panels.";
        evidence = "Discussion: 'Under the scenario of a 10-fold increased k cleave in the tumor relative to the periphery, Figure 4d suggests ...'; the supplement snapshot also has 10x. Whether the same fold applies to Figure 4a-c is not stated." * PRE,
        basis = :interpretation, stage = :before_comparison, figures = ["fig4"]),
    Fix("breathing_kopen_supplement", :ambiguity,
        "Mask breathing: kopen = 1.16e-2 1/s and kclose = Kmask x kopen (the main text gives only Kmask = kclose/kopen).";
        evidence = "Supplement s003 parameter list: kopen 1.16e-2, kclose 0.232 = 20 x 1.16e-2 for Kmask 20. The other split (kclose fixed) is a sensitivity run and is immaterial (report.md)." * PRE,
        basis = :supplement, stage = :before_comparison),
    Fix("dose_as_intact", :ambiguity,
        "The Pb-Tx dose enters plasma as the intact form with both masks closed (c_c); the parental mAb0,0 as m_m.";
        evidence = "Results: 'The Pb-Tx is depicted first as an intact moiety with two masks'. Whether the masks start closed or at the breathing equilibrium is not stated (breathing equilibrates within seconds, so this has no visible effect)." * PRE,
        basis = :interpretation, stage = :before_comparison),
    Fix("receptor_initial", :missing,
        "Peripheral and tumour receptor start at steady state ksynR/(kendo V); the central receptor state starts at 0 with no synthesis.";
        evidence = "The initial receptor levels are not stated; s003 lists a central receptor state but prints no zero-order (synthesis) vector. The central choice is in notes.md before comparing; the steady-state start is not recorded there, so the stage is not known and is set to :after_comparison.",
        basis = :assumption, stage = :after_comparison),
    Fix("readout_total_all_forms", :ambiguity,
        "'Total Drug' (nM) is the sum of all unbound and receptor-bound forms in the compartment (excluding internalised and eliminated pools). In the fixed run, Figure 4b,c use `unbound_readout` instead.";
        evidence = "Figure 3 and Figure 4a-c plot 'Total Drug (nM)'; the paper does not define which forms it sums. In plasma both readings are identical (no central receptor)." * PRE,
        basis = :interpretation, stage = :before_comparison),
    Fix("readout_uptake_internalised", :ambiguity,
        "Figure 4d,e 'Total Drug (nmol)' is the cumulative receptor-mediated uptake: the internalised pools tmdd_m + tmdd_o times the compartment volume.";
        evidence = "Figure 4d,e are described in the Discussion as receptor-mediated uptake in periphery and tumour; the y label says 'Total Drug (nmol)'." * PRE,
        basis = :interpretation, stage = :before_comparison, figures = ["fig4"]),
]

const FIXES = [
    PUBLISHED_FILLS...,
    Fix("unbound_readout", :from_target,
        "Figure 4b,c periphery/tumour 'Total Drug' is read as the sum of the six UNBOUND forms. The as-published run used the other reading (unbound + target-bound, the same sum as for plasma).";
        evidence = "The paper does not define which forms 'Total Drug' sums in periphery and tumour (neither reading is the published one). The unbound reading was chosen because it matches Fig 4b,c for mAb0,0 and 10x (e.g. mAb ~1e-2 nM in tumour, far below the 270 nM target level, Table 2); selected by comparison with the figure, hence circular.",
        basis = :figure, stage = :after_comparison, figures = ["fig4"]),
]

const MONKEY_BW = 2.6      # kg, Table 1
const HUMAN_BW = 70.0      # kg, Table 2
const HUMAN_DOSE = 4.5     # mg/kg, Figure 4 caption
const MASKS = Dict("M1" => 220.0, "M2" => 57.0)             # Table 1 fold-masking Kmask
const KCLEAVE_MAX = 3.0e-7                                   # Table 1/2 '<3e-7': the bound, used as the value (fixes kcleave_monkey_3e-7, kcleave_human_3e-7)
const HUMAN_MASKS = Dict("mAb" => nothing, "10x" => 10.0, "100x" => 100.0, "1000x" => 1000.0, "10000x" => 10000.0)

"""
Simulation specs. `key` identifies one ODE run; `out` is the observed variable the figure plots.
Fig 3: plasma total drug for a molecule at a dose. Fig 4: one run per mask strength, five outputs.
"""
struct SimSpec
    species::String       # "monkey" or "human"
    molecule::String      # "mAb", "M1S1", ..., "10x", ...
    dose::Float64         # mg/kg
    tend::Float64         # days
end

const FIG3 = Dict(
    "fig3a" => (3.0, ["mAb00", "M1S1", "M2S2", "M1_0"]),
    "fig3b" => (5.0, ["mAb00", "M2S1", "M1S2"]),
    "fig3c" => (10.0, ["mAb00", "M2S1", "M1S2"]))
const FIG4_OUT = Dict("a" => :Ctot_central, "b" => :Ctot_peripheral, "c" => :Ctot_tumor, "d" => :Upt_peripheral, "e" => :Upt_tumor)

function spec_for(id::AbstractString)
    m = match(r"^fig3([abc])_(.+)$", id)
    if m !== nothing
        return SimSpec("monkey", m[2], FIG3["fig3" * m[1]][1], 28.0), :Ctot_central
    end
    m = match(r"^fig4([a-e])_(.+)$", id)
    m === nothing && error("unknown reference id $id")
    return SimSpec("human", m[2], HUMAN_DOSE, 56.0), FIG4_OUT[m[1]]
end

"Parameters, dosed state and mask for a molecule. Each unstated input is taken from its fix (NotRunnable without it)."
function molecule_setup(sp::SimSpec; fixes = Fix[])
    has(id) = any(f -> f.id == id, fixes)
    need(id, what) = has(id) || throw(NotRunnable(what))
    need("dose_as_intact", "the dosed molecular form is not stated")
    kclose_fixed = has("breathing_kclose_fixed")
    kclose_fixed || need("breathing_kopen_supplement", "the absolute breathing rates kopen / kclose are not stated (only Kmask)")
    if sp.species == "monkey"
        if sp.molecule == "mAb00"
            Km, kcl, dosed = 1.0, 0.0, :pb_m_m_central
        else
            mask, sub = sp.molecule[1:2], sp.molecule[3:end]
            Km = MASKS[mask]; dosed = :pb_c_c_central
            kcl = sub == "_0" ? 0.0 : (need("kcleave_monkey_3e-7", "monkey kcleave is given only as '<3e-7'"); KCLEAVE_MAX)
        end
        kopen = kclose_fixed ? 0.232 / Km : KOPEN
        return monkey_pars(; Kmask = Km, kcleave = kcl, kopen), dosed, MONKEY_BW
    else
        if sp.molecule == "mAb"
            Km, dosed = 1.0, :pb_m_m_central
        else
            Km, dosed = HUMAN_MASKS[sp.molecule], :pb_c_c_central
        end
        need("kcleave_human_3e-7", "human kcleave is given only as '<3e-7'")
        need("tumour_kcleave_10fold", "the tumour kcleave fold is not stated for all of Figure 4")
        need("human_V2_2.6L", "human V2 is not in Table 2")
        kopen = kclose_fixed ? 0.232 / Km : KOPEN
        return human_pars(; Kmask = Km, kopen, kcleave = KCLEAVE_MAX, kcleave_tumor_fold = 10.0, V2 = 2.6), dosed, HUMAN_BW
    end
end

const _SYS = Ref{Any}(nothing)
get_system() = (_SYS[] === nothing && (_SYS[] = build_system()); _SYS[])

"ODE problem plus dosing callback for one spec. Returns (prob, callback, saveat)."
function make_problem(sp::SimSpec; fixes = Fix[], npts = 1500)
    sys = get_system()
    pars, dosed, bw = molecule_setup(sp; fixes)
    any(f -> f.id == "molecular_weight_150kDa", fixes) || throw(NotRunnable("the molecular weight for the mg/kg to nM conversion is not stated"))
    any(f -> f.id == "receptor_initial", fixes) || throw(NotRunnable("the initial receptor levels are not stated"))
    amount = dose_nM(sp.dose, bw, pars[:V1])
    prob = ODEProblem(sys, problem_pairs(sys, pars), (0.0, sp.tend))
    cb = bolus_callback(sys, getproperty(sys, dosed), [0.0], amount)
    return prob, cb, collect(range(0.0, sp.tend; length = npts))
end

const _CACHE = Dict{Any, Any}()
function run_spec(sp::SimSpec; fixes = Fix[])
    key = (sp, [f.id for f in fixes])
    get!(_CACHE, key) do
        prob, cb, ts = make_problem(sp; fixes)
        sol = solve(prob, Rodas5P(); callback = cb, saveat = ts, abstol = 1e-10, reltol = 1e-8)
        sol.retcode == ReturnCode.Success || error("solver failed for $(sp): $(sol.retcode)")
        (sys = get_system(), sol = sol, ts = ts)
    end
end

"`(x [day], y)` for reference curve `ref`."
function simulate_curve(ref; fixes = Fix[])
    sp, out = spec_for(ref.id)
    has(id) = any(f -> f.id == id, fixes)
    if startswith(String(out), "Upt_")
        has("readout_uptake_internalised") || throw(NotRunnable("Figure 4d,e 'Total Drug (nmol)' is not defined"))
    elseif has("unbound_readout") && out in (:Ctot_peripheral, :Ctot_tumor)
        out = Symbol(replace(String(out), "Ctot_" => "Cfree_"))
    else
        has("readout_total_all_forms") || throw(NotRunnable("'Total Drug' is not defined"))
    end
    r = run_spec(sp; fixes)
    return r.ts, r.sol(r.ts; idxs = getproperty(r.sys, out)).u
end

"Problems for the solver benchmarks: monkey 3 mg/kg Pb-Tx M1,S1 and human 4.5 mg/kg 1000x mask (stiff, 55 states)."
function benchmark_setups()
    out = BenchmarkSetup[]
    for (name, sp) in [("monkey_3mgkg_M1S1", SimSpec("monkey", "M1S1", 3.0, 28.0)), ("human_4.5mgkg_1000x", SimSpec("human", "1000x", HUMAN_DOSE, 56.0))]
        prob, cb, ts = make_problem(sp; fixes = FIXES)
        push!(out, BenchmarkSetup(prob; name, callback = cb, saveat = ts))
    end
    return out
end
