# Protocol, targets, tolerances, fixes for ouerdani_2015 (Figure 4, typical simulations).
# TARGETS FROZEN 2026-11-04 09:44 IST: reference/*.csv, digitise.toml and the tolerances below
# were written before any simulation was compared with them (no simulation had been run).
using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq

# Default tolerance (RMS 0.015, P95 0.03 of the axis box), no per-curve change.
const TOLERANCE = Tolerance()
const TOLERANCES = Dict{String, Tolerance}()

# Days per month used to put the clinical model time (days) on Figure 4's month axis.
# Not stated in the paper; 365.25/12 assumed, 30 and 30.4375 differ by 1.5 %.
const DAYS_PER_MONTH = 365.25 / 12

# The "as published" choices below are the readings available from the text alone; the
# paper does not state P0 or which dose the mouse panel uses (see FIXES and report.md).
#   mice:     P0 = 175 mm^3 (midpoint of the 100-250 mm^3 inclusion range, Methods),
#             AUC = 656.8 (30 mg/kg, middle dose), horizon 75 d (Figure 4 axis).
#   patients: P0 = 75 mm (rough median baseline SLD in Figure 1 right; the paper gives none),
#             AUC = 771.6 ug*h/mL (stated mean exposure), horizon 82 months (axis).
const TEXT_ONLY = (mouse = (P0 = 175.0, AUC = MOUSE_AUC.dose30), patient = (P0 = 75.0, AUC = PATIENT_AUC))

# Values used by the fixes (found by trying the readings, see report.md):
#   p0: P0 is not given anywhere ("set to the observed value", per individual). Taken from
#       Figure 4 at t = 0 (150 mm^3; 66.8 mm), which is circular but the only source.
#   mouse_auc: the typical mouse panel does not say which dose group it shows; the
#       100 mg/kg exposure (AUC 1140.8) reproduces it, 10 and 30 mg/kg do not.
const FIXED_P0 = (mouse = 150.0, patient = 66.8)
const FIXED_MOUSE_AUC = MOUSE_AUC.dose100

const FIXES = [
    Fix("p0", :from_target,
        "Initial tumour size P0 is not stated for the typical simulation; set to 150 mm^3 (mice) and 66.8 mm (patients), read from Figure 4 at t = 0.";
        evidence = "Methods: P0 'set to the observed value' (individual baselines); no typical value given in text, Table 1, s001.txt/s004.txt (P0 = BASE/VB data column). Value read from the curve start in PSP4-4-660-g004.jpg (circular)."),
    Fix("mouse_auc", :from_target,
        "Mouse panel simulated with the 100 mg/kg exposure (AUC 1140.8 ug*h/mL); the paper does not say which dose group the typical mouse simulation uses.";
        evidence = "Only gamma depends on AUC (beta_alpha fixed to 0). With P0 = 150: 100 mg/kg rms 0.001, 30 mg/kg rms 0.074, 10 mg/kg rms 0.27 (explore log)."),
]

function setting(ref_id::AbstractString, fixes)
    ids = Set(f.id for f in fixes)
    mouse = startswith(ref_id, "fig4_left")
    d = mouse ? TEXT_ONLY.mouse : TEXT_ONLY.patient
    P0 = ("p0" in ids) ? (mouse ? FIXED_P0.mouse : FIXED_P0.patient) : d.P0
    AUC = (mouse && "mouse_auc" in ids) ? FIXED_MOUSE_AUC : d.AUC
    return (mouse = mouse, P0 = P0, AUC = AUC)
end

function problem_for(ref_id; fixes = Fix[], P0 = nothing, AUC = nothing, par = nothing)
    s = setting(ref_id, fixes)
    P0 = P0 === nothing ? s.P0 : P0
    AUC = AUC === nothing ? s.AUC : AUC
    par = par === nothing ? (s.mouse ? PRECLINICAL : CLINICAL) : par
    tspan = s.mouse ? (0.0, 75.0) : (0.0, 82.0 * DAYS_PER_MONTH)
    return tumour_problem(par, P0, AUC, tspan), s.mouse
end

function simulate_curve(ref; fixes = Fix[], kw...)
    prob, mouse = problem_for(ref.id; fixes, kw...)
    tend = prob.tspan[2]
    ts = collect(range(0.0, tend; length = 1000))
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10, saveat = ts)
    state = endswith(ref.id, "_P") ? sol[TUMOUR.P] : sol[TUMOUR.K]
    return (mouse ? ts : ts ./ DAYS_PER_MONTH), state
end

function benchmark_setups()
    out = BenchmarkSetup[]
    for (name, key) in (("mouse typical, Figure 4 left", "fig4_left_P"), ("patient typical, Figure 4 right", "fig4_right_P"))
        prob, mouse = problem_for(key; fixes = FIXES)
        ts = collect(range(0.0, prob.tspan[2]; length = 1000))
        push!(out, BenchmarkSetup(; name, prob, saveat = ts))
    end
    return out
end
