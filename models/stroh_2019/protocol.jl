# Protocol for stroh_2019: dosing, simulation settings, targets, tolerances, fixes, benchmark setups.
using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq

# ---- targets (frozen 2026-11-04 09:48 local, before any simulation was compared; see notes.md)
# Grade-bearing: Figure 3a-c (monkey calibration lines, 10 curves) and Figure 4a-e (human projections, 25 curves).
# Not graded: observed monkey points, Figure 5.
const TOLERANCE = Tolerance()                       # default rms 0.015, p95 0.03 of the axis box
const TOLERANCES = Dict{String, Tolerance}()        # no per-curve overrides

# ---- fixes. The two `:missing` gap-fills (molecular weight, human V2) are needed for the model to run at all, so the
# "as published" run uses them too (same values in both runs); they are listed here so that the report states how each gap is
# filled. Only `unbound_readout` changes the simulated curves between the two runs.
const FIXES = [
    Fix("molecular_weight_150kDa", :missing,
        "mg/kg doses are converted to nM with MW = 150 kDa (model.jl MW). Applied in BOTH the as-published and the fixed run.";
        evidence = "The paper gives no molecular weight. 150 kDa is the standard IgG1 mass (Pb-Tx adds a prodomain, so the true mass is somewhat higher; value is a standard constant, not read from the figure). Consistency only, not the source: 4.5 mg/kg x 70 kg / 150 kDa / 2.6 L = 808 nM."),
    Fix("human_V2_2.6L", :missing,
        "Human peripheral volume V2 = 2.6 L (not in Table 2). Applied in BOTH the as-published and the fixed run.";
        evidence = "Derived from the paper's own tables: Table 1 monkey receptor concentration 9.0e-5 nmol/s / (1e-4 1/s * 0.1 L) = 9 nM; Results state normal-tissue receptor concentration is assumed equal in human; Table 2 ksynR peripheral 2.3e-3 nmol/s / (1e-4 * 9 nM) = 2.56 L ~ 2.6 L (= V1)."),
    Fix("unbound_readout", :from_target,
        "Figure 4b,c periphery/tumour 'Total Drug' is read as the sum of the six UNBOUND forms. The as-published run used the other reading (unbound + target-bound, the same sum as for plasma).";
        evidence = "The paper does not define which forms 'Total Drug' sums in periphery and tumour (neither reading is the published one). The unbound reading was chosen because it matches Fig 4b,c for mAb0,0 and 10x (e.g. mAb ~1e-2 nM in tumour, far below the 270 nM target level, Table 2); selected by comparison with the figure, hence circular."),
]

const MONKEY_BW = 2.6      # kg, Table 1
const HUMAN_BW = 70.0      # kg, Table 2
const HUMAN_DOSE = 4.5     # mg/kg, Figure 4 caption
const MASKS = Dict("M1" => 220.0, "M2" => 57.0)             # Table 1 fold-masking Kmask
const KCLEAVE_MAX = 3.0e-7                                   # Table 1: upper bound carried forward for S1, S2
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

"Parameters, dosed state and mask for a molecule. Fix id `kclose_fixed` switches the breathing split (see report)."
function molecule_setup(sp::SimSpec; fixes = Fix[])
    kclose_fixed = any(f -> f.id == "breathing_kclose_fixed", fixes)
    if sp.species == "monkey"
        if sp.molecule == "mAb00"
            Km, kcl, dosed = 1.0, 0.0, :pb_m_m_central
        else
            mask, sub = sp.molecule[1:2], sp.molecule[3:end]
            Km = MASKS[mask]; kcl = sub == "_0" ? 0.0 : KCLEAVE_MAX; dosed = :pb_c_c_central
        end
        kopen = kclose_fixed ? 0.232 / Km : KOPEN
        return monkey_pars(; Kmask = Km, kcleave = kcl, kopen), dosed, MONKEY_BW
    else
        if sp.molecule == "mAb"
            Km, dosed = 1.0, :pb_m_m_central
        else
            Km, dosed = HUMAN_MASKS[sp.molecule], :pb_c_c_central
        end
        kopen = kclose_fixed ? 0.232 / Km : KOPEN
        return human_pars(; Kmask = Km, kopen), dosed, HUMAN_BW
    end
end

const _SYS = Ref{Any}(nothing)
get_system() = (_SYS[] === nothing && (_SYS[] = build_system()); _SYS[])

"ODE problem plus dosing callback for one spec. Returns (prob, callback, saveat)."
function make_problem(sp::SimSpec; fixes = Fix[], npts = 1500)
    sys = get_system()
    pars, dosed, bw = molecule_setup(sp; fixes)
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
    if any(f -> f.id == "unbound_readout", fixes)
        out = Symbol(replace(String(out), "Ctot_" => "Cfree_"))
    end
    r = run_spec(sp; fixes)
    return r.ts, r.sol(r.ts; idxs = getproperty(r.sys, out)).u
end

"Problems for the solver benchmarks: monkey 3 mg/kg Pb-Tx M1,S1 and human 4.5 mg/kg 1000x mask (stiff, 55 states)."
function benchmark_setups()
    out = BenchmarkSetup[]
    for (name, sp) in [("monkey_3mgkg_M1S1", SimSpec("monkey", "M1S1", 3.0, 28.0)), ("human_4.5mgkg_1000x", SimSpec("human", "1000x", HUMAN_DOSE, 56.0))]
        prob, cb, ts = make_problem(sp)
        push!(out, BenchmarkSetup(prob; name, callback = cb, saveat = ts))
    end
    return out
end
