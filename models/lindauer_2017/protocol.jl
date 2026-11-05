# Protocol for the lindauer_2017 reproduction: Figure 4 (deterministic one-at-a-time sensitivity analysis of the
# human melanoma model, 26 weeks, 5 and 0.5 mg/kg). Targets and tolerances were frozen on 2026-11-04 09:48 IST
# before any simulation was run (see notes.md).
using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

# Default tolerance (rms 0.015, p95 0.03 of the axis box), no per-curve overrides.
const TOLERANCE = Tolerance()
const TOLERANCES = Dict{String, Tolerance}()

# --- choices the paper does not state (as-published values declared in the notes before running) -------------
const BW_KG = 70.0            # human body weight used for the allometric scaling in the Methods
const DOSE_INTERVAL_DAYS = 21.0   # Q3W (Figure 3 is Q3W; Figure 4 schedule is not stated)
const T_END_DAYS = 182.0      # 26 weeks
const DOSE_MGKG = Dict("5mgkg" => 5.0, "0p5mgkg" => 0.5)

"Settings of one run, built from its fixes by `setup_from_fixes`."
struct Setup
    scenario::String
    V_blood::Float64
    interval_days::Float64
    KdegPD1::Float64
end
Setup(; scenario = "fast_allom", V_blood = 1400.0, interval_days = DOSE_INTERVAL_DAYS, KdegPD1 = 0.00246) =
    Setup(scenario, V_blood, interval_days, KdegPD1)

# Provenance ledger (2026-11-05): every choice the paper does not state is a Fix in the run that uses it. The fills
# below were declared in notes.md before any simulation ("As-published run (declared a priori)"). Without the scenario,
# V_blood, schedule, body-weight and bolus fills the simulation throws NotRunnable.
const PRE = " Declared in notes.md before any simulation was run."
const SHARED_FILLS = [
    Fix("schedule_Q3W", :missing,
        "Figure 4 doses every 3 weeks (bolus at days 0, 21, ..., 168; response at 26 weeks).";
        evidence = "The Figure 4 schedule is not stated. Q3W is the schedule of Figure 3a ('once every 3 weeks') and the regimen the paper recommends (2 mg/kg Q3W)." * PRE,
        basis = :interpretation, stage = :before_comparison),
    Fix("dose_body_weight_70kg", :missing,
        "mg/kg doses converted to amounts with a 70 kg body weight (MW 149000 g/mol from s006).";
        evidence = "Methods: 'A standard body weight of 20 g and 70 kg was assumed for mice and humans' (stated for the allometric scaling; its use for dosing is not stated)." * PRE,
        basis = :interpretation, stage = :before_comparison),
    Fix("bolus_doses", :missing,
        "Each human dose is an IV bolus into the central compartment.";
        evidence = "Dose administration in the human simulations (bolus or infusion, duration) is not stated; no effect at the 26-week readout (report.md, Issues)." * PRE,
        basis = :assumption, stage = :before_comparison),
    Fix("table_s2_full_precision", :ambiguity,
        "Species-invariant estimates taken at Table S2 precision (gamma 2.28, Tmulti 4.32, EMAXTP 94.7, EC50TP 1.46) instead of the Table 1 roundings (gamma 2.3, Tmulti 4.3).";
        evidence = "Table 1 and Table S2 print the same estimates at different precision; Table S2 is the estimation output. Table 1 roundings give a base case of -1.340 / -0.418 instead of -1.302 / -0.401 (notes.md)." * PRE,
        basis = :supplement, stage = :before_comparison),
    Fix("equations_from_s006", :code,
        "Where the supplement equations (s001) and the NONMEM control stream (s006) differ, s006 is followed: Ceb in both endosomal transfer terms, the full central-compartment equation, K_IVIV = 1, MW 149000. The s006 branch DE = SLtg*RO for RO <= 0.01 % is not applied (no effect, report.md).";
        evidence = "s001 prints dCeb/dt and dFcRn/dt with Cvs where s006 and mass balance have Ceb, and a truncated central-compartment equation (report.md, Issues)." * PRE,
        basis = :supplement, stage = :before_comparison),
]

const PUBLISHED_FILLS = [
    Fix("scenario_table_s3_row1", :ambiguity,
        "Figure 4 base case = Table S3 row 1, fast growth with allometric SLtg (L0 0.0088/day, SLtg 2.575E-06).";
        evidence = "The Fig S1 flowchart (s003) says the sensitivity analysis takes 'one scenario (melanoma, 1)'; Table S3 lists fast / allometric first. The Fig 3a legend order would make row 4 'scenario 1' instead (report.md, Issues)." * PRE,
        basis = :interpretation, stage = :before_comparison),
    Fix("human_V_blood_mouse_value", :missing,
        "Human V_blood (blood volume entering the PD-1 target amount in blood) = 1400 uL, the only V_blood printed (the mouse value in s006).";
        evidence = "Table 1 has no V_blood row and no human value is printed; s006 'TVV_blood = THETA(18)' = 1400 uL for a 20 g mouse. Nothing says it is used for man." * PRE,
        basis = :assumption, stage = :before_comparison),
    SHARED_FILLS...,
]

"""
Fixes of the "with fixes" run: three changes against the as-published run plus the shared fills. The candidate menu
(V_blood, scenario row, schedule) was fixed in the notes before any simulation, and the scenario was chosen on the two
caption base values only (identification grid in notes.md).
"""
const FIXES = [
    Fix("human_V_blood_5L", :missing,
        "Human V_blood (blood volume entering the PD-1 target amount in blood) set to 5 L = 5e6 uL. The paper gives no human value and Table 1 has no V_blood row; 1400 uL is the mouse blood volume (s006 THETA(18)).";
        evidence = "Standard adult blood volume (70 kg x 70 mL/kg). s006 line 'TVV_blood = THETA(18) ; Blood volume' = 1400 (20 g mouse); Figure 4/Results list V_blood among the sensitive human parameters; with 1400 uL the base case is -6.0 / -4.7 (log ratio) instead of -1.29 / -0.39. " *
                   "The 5 L value was declared in notes.md before any simulation (candidate menu), not tuned (5.2 L would fit the caption slightly better).",
        basis = :standard, stage = :before_comparison),
    Fix("KdegPD1_allometric", :typo,
        "Human KdegPD1 = 0.0194 * (70 kg / 0.02 kg)^-0.25 = 0.002522 1/h instead of the printed Table 1 value 0.00246. The Methods state this formula; the printed value is 2.5 % lower and cannot be obtained from it (SLtg, scaled by the same rule, is reproduced exactly: 2.575E-06).";
        evidence = "Methods, equation for KdegPD1_human (BW 20 g / 70 kg, exponent -0.25) vs Table 1 'Value in man' 0.00246; Table S3 SLtg 2.575E-06 = 1.98E-05 * 3500^-0.25. " *
                   "Noted as an inconsistency before running, applied after seeing its effect on the base case; the evidence is textual and independent of Figure 4.",
        basis = :paper, stage = :after_comparison),
    Fix("scenario_med_growth", :from_target,
        "Base case of Figure 4 = Table S3 row 4: medium growth rate (L0 0.0036/day) with growth-proportional scaling of SLtg (6.247E-07). The paper does not say which of the six scenarios is the base case (Fig S1: 'one scenario (melanoma, 1)').";
        evidence = "Of the 12 scenario x schedule combinations only medium growth / growth-proportional, Q3W, gives the Figure 4 caption base values (-1.302 / -0.401 vs printed -1.29 / -0.39); next best -1.44 / -0.46. Chosen by matching fig4_base_case.",
        basis = :figure, stage = :after_comparison),
    SHARED_FILLS...,
]

function setup_from_fixes(fixes)
    ids = Set(f.id for f in fixes)
    scenario = "scenario_med_growth" in ids ? "med_growth" : "scenario_table_s3_row1" in ids ? "fast_allom" :
        throw(NotRunnable("the Figure 4 base scenario (row of Table S3) is not stated"))
    V_blood = "human_V_blood_5L" in ids ? 5.0e6 : "human_V_blood_mouse_value" in ids ? 1400.0 :
        throw(NotRunnable("no human V_blood is printed"))
    "schedule_Q3W" in ids || throw(NotRunnable("the Figure 4 dosing schedule is not stated"))
    ("dose_body_weight_70kg" in ids && "bolus_doses" in ids) || throw(NotRunnable("dose amount and route of the human simulations are not stated"))
    KdegPD1 = "KdegPD1_allometric" in ids ? 0.0194 * (BW_KG / 0.02)^(-0.25) : 0.00246
    return Setup(; scenario, V_blood, interval_days = DOSE_INTERVAL_DAYS, KdegPD1)
end

# --- tornado parameters: figure label -> raw parameter field -----------------------------------------------
const TORNADO_FIELD = Dict("gamma" => :gamma, "FR" => :FR, "SLtg" => :SLtg, "L0" => :L0, "Tmulti" => :Tmulti,
    "N_PD1_TC" => :N_PD1_TC, "V_blood" => :V_blood, "N_Tcell" => :N_Tcell, "EMAXTP" => :EMAXTP, "V1" => :V1,
    "KdegPD1" => :KdegPD1, "CL" => :CL, "PLQ" => :PLQ, "L" => :Lfrac, "CLup" => :CLup, "V2" => :V2, "EC50TP" => :EC50TP)

const _SYS = Ref{Any}(nothing)
lindauer_system() = (_SYS[] === nothing && (_SYS[] = build_lindauer()); _SYS[])

dose_times_h(interval_days) = repeated_doses(0.0, interval_days * 24, ceil(Int, T_END_DAYS / interval_days))
dose_nmol(mgkg) = mgkg * BW_KG * 1e6 / MW   # mg -> ug (x1000) -> nmol (x1000/MW)

"ODEProblem + dosing callback for one human run."
function human_problem(raw, mgkg, interval_days)
    sys = lindauer_system()
    times = dose_times_h(interval_days)
    cb = bolus_callback(sys, sys.A1, times, dose_nmol(mgkg))
    prob = ODEProblem(sys, problem_pairs(sys, raw), (0.0, T_END_DAYS * 24); callback = cb)
    return prob, cb
end

const SOLVE_KW = (abstol = 1e-14, reltol = 1e-9)
function solve_human(raw, mgkg, interval_days; saveat = nothing, alg = Rodas5P())
    prob, cb = human_problem(raw, mgkg, interval_days)
    sys = lindauer_system()
    kw = saveat === nothing ? (; save_everystep = false, save_start = false) : (; saveat)
    sol = solve(prob, alg; SOLVE_KW..., kw..., dtmax = 24.0)
    return sol, sys
end

"log10(tumour volume at 26 weeks / initial tumour volume)."
function log_ratio(raw, mgkg, interval_days)
    sol, sys = solve_human(raw, mgkg, interval_days)
    # a failed solve (retcode not successful) is NaN, never a finite stand-in value; the scorer counts NaN as the maximum error
    SciMLBase.successful_retcode(sol) || return NaN
    v = sol(T_END_DAYS * 24; idxs = sys.TV)
    return v > 0 ? log10(v / raw.W0) : NaN
end

const _CACHE = Dict{Any, Float64}()
function cached_log_ratio(s::Setup, mgkg, name, mult)
    return get!(_CACHE, (s, mgkg, name, mult)) do
        raw = human_raw(; scenario = s.scenario, V_blood = s.V_blood, KdegPD1 = s.KdegPD1)
        if name != "base"
            f = TORNADO_FIELD[name]
            raw = merge(raw, NamedTuple{(f,)}((getfield(raw, f) * mult,)))
        end
        log_ratio(raw, mgkg, s.interval_days)
    end
end

"Simulated bar values for one reference id (fig4_<dose>_x5, fig4_<dose>_x0p2)."
function simulate_bars(ref, fixes)
    s = setup_from_fixes(fixes)
    m = match(r"^fig4_(5mgkg|0p5mgkg)_(x5|x0p2)$", ref.id)
    m === nothing && error("unknown reference id $(ref.id)")
    order = String.(split(ref.meta["order"], ", "))
    mult = m[2] == "x5" ? 5.0 : 0.2
    return [cached_log_ratio(s, DOSE_MGKG[m[1]], name, mult) for name in order]
end

"Simulated values at exactly the reference x positions (kind: bars scoring)."
function simulate_curve(ref; fixes = Fix[])
    if ref.id == "fig4_base_case"
        s = setup_from_fixes(fixes)
        return ref.x, [cached_log_ratio(s, DOSE_MGKG["5mgkg"], "base", 1.0), cached_log_ratio(s, DOSE_MGKG["0p5mgkg"], "base", 1.0)]
    end
    return ref.x, simulate_bars(ref, fixes)
end

"""
Exact vertical errors of the bar encoding (simulated minus digitised, log10 units and fraction of the axis box).
Same quantity as the scorer's `:bars` comparison, plus signed values (NaN for failed solves).
"""
function vertical_errors(ref; fixes = Fix[])
    ys = ref.id == "fig4_base_case" ? [cached_log_ratio(setup_from_fixes(fixes), DOSE_MGKG["5mgkg"], "base", 1.0),
        cached_log_ratio(setup_from_fixes(fixes), DOSE_MGKG["0p5mgkg"], "base", 1.0)] : simulate_bars(ref, fixes)
    box = ref.yaxis.hi - ref.yaxis.lo
    e = ys .- ref.y
    return (; err = e, frac = e ./ box, rms = sqrt(sum(abs2, e ./ box) / length(e)), maxabs = maximum(abs, e ./ box))  # NaN propagates for failed solves
end

"ODE problems used for the figures, for the solver benchmarks (base case, as published, 5 and 0.5 mg/kg Q3W, 26 weeks)."
function benchmark_setups()
    out = BenchmarkSetup[]
    for (nm, mg) in (("5mgkg", 5.0), ("0p5mgkg", 0.5))
        prob, cb = human_problem(human_raw(), mg, DOSE_INTERVAL_DAYS)
        push!(out, BenchmarkSetup(prob; name = "lindauer_2017 human base, $(mg) mg/kg Q3W, 26 weeks", callback = cb,
            saveat = collect(range(0.0, T_END_DAYS * 24; length = 366))))
    end
    return out
end
