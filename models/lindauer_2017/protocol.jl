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

"Settings of the as-published run; fixes override single entries."
struct Setup
    scenario::String
    V_blood::Float64
    interval_days::Float64
    KdegPD1::Float64
end
Setup(; scenario = "fast_allom", V_blood = 1400.0, interval_days = DOSE_INTERVAL_DAYS, KdegPD1 = 0.00246) =
    Setup(scenario, V_blood, interval_days, KdegPD1)

"""
Fixes. Both change a protocol value that the paper does not state; the candidate menu (V_blood, scenario row, schedule) was
fixed in the notes before any simulation, and the choice was made on the two caption base values only (identification grid in
notes.md). The 64 bars of Figure 4 are then an independent check.
"""
const FIXES = [
    Fix("human_V_blood_5L", :missing,
        "Human V_blood (blood volume entering the PD-1 target amount in blood) set to 5 L = 5e6 uL. The paper gives no human value and Table 1 has no V_blood row; 1400 uL is the mouse blood volume (s006 THETA(18)).";
        evidence = "s006 line 'TVV_blood = THETA(18) ; Blood volume' = 1400 (20 g mouse); Figure 4/Results list V_blood among the sensitive human parameters; with 1400 uL the base case is -6.0 / -4.7 (log ratio) instead of -1.29 / -0.39 (identification grid in notes.md)."),
    Fix("KdegPD1_allometric", :typo,
        "Human KdegPD1 = 0.0194 * (70 kg / 0.02 kg)^-0.25 = 0.002522 1/h instead of the printed Table 1 value 0.00246. The Methods state this formula; the printed value is 2.5 % lower and cannot be obtained from it (SLtg, scaled by the same rule, is reproduced exactly: 2.575E-06).";
        evidence = "Methods, equation for KdegPD1_human (BW 20 g / 70 kg, exponent -0.25) vs Table 1 'Value in man' 0.00246; Table S3 SLtg 2.575E-06 = 1.98E-05 * 3500^-0.25."),
    Fix("scenario_med_growth", :from_target,
        "Base case of Figure 4 = Table S3 row 4: medium growth rate (L0 0.0036/day) with growth-proportional scaling of SLtg (6.247E-07). The paper does not say which of the six scenarios is the base case (Fig S1: 'one scenario (melanoma, 1)').";
        evidence = "Of the 12 scenario x schedule combinations only medium growth / growth-proportional, Q3W, gives the Figure 4 caption base values (-1.302 / -0.401 vs printed -1.29 / -0.39); next best -1.44 / -0.46."),
]

function setup_from_fixes(fixes)
    s = Dict{Symbol, Any}()
    for f in fixes
        if f.id == "human_V_blood_5L"
            s[:V_blood] = 5.0e6
        elseif startswith(f.id, "scenario_")
            s[:scenario] = f.id[length("scenario_")+1:end]
        elseif f.id == "schedule_Q2W"
            s[:interval_days] = 14.0
        elseif f.id == "KdegPD1_allometric"
            s[:KdegPD1] = 0.0194 * (BW_KG / 0.02)^(-0.25)
        end
    end
    return Setup(; s...)
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
