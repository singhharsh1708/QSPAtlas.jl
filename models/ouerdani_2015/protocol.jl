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

# Every input the paper does not state is a Fix, in both runs (provenance ledger, 2026-11-05).
# The "as published" run needs four fills to run at all; without them `simulate_curve` throws NotRunnable.
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

const STAGE_NOTE = " Stage: the notes record only that these inputs were tried after the 09:44 freeze; the exact order is not recorded, so :after_comparison."

# fills shared by both runs
const PATIENT_AUC_FIX = Fix("patient_auc_mean", :ambiguity,
    "Patient panel simulated with the stated mean exposure AUC = 771.6 ug*h/mL; the caption does not say which exposure Figure 4 right uses.";
    evidence = "Results: mean exposure of the 47 patients 771.6 ug*h/mL (range 629.4-802.4); Supplementary Table 2 uses 800. 771.6 is the only value the text ties to the population. Not chosen by matching: 800 fits K about six times tighter and was not adopted (report.md, Sensitivity)." * STAGE_NOTE,
    basis = :interpretation, stage = :after_comparison, figures = ["fig4"])
const DAYS_PER_MONTH_FIX = Fix("days_per_month", :missing,
    "Clinical model time (days) put on the Figure 4 month axis with 365.25/12 = 30.4375 days per month.";
    evidence = "Table 1 and the NONMEM code use days; Figure 4 right is in months and the conversion is not stated. 365.25/12 is the mean calendar month (30 also passes, 28 fails; report.md)." * STAGE_NOTE,
    basis = :standard, stage = :after_comparison, figures = ["fig4"])
const CONTINUOUS_FIX = Fix("continuous_treatment", :ambiguity,
    "Treatment is continuous over the whole horizon (75 d mice, 82 months patients): constant AUC and the exp(-delta t) term from t = 0, no stop.";
    evidence = "Eq. 1-3 and the NONMEM \$DES blocks (s001, s004) have no dosing stop; the caption does not say how long treatment lasts (the mice were dosed for 24 days only)." * STAGE_NOTE,
    basis = :interpretation, stage = :after_comparison, figures = ["fig4"])

# The as-published run: the readings available from the text alone, plus the shared fills.
const PUBLISHED_FILLS = [
    Fix("p0_text_guess", :missing,
        "Initial tumour size P0 for the typical simulation: 175 mm^3 (mice, midpoint of the 100-250 mm^3 inclusion range) and 75 mm (patients, rough median baseline SLD read from Figure 1 right).";
        evidence = "The paper gives no typical P0 (Methods: 'set to the observed value'). 175 is the midpoint of the Methods inclusion range; 75 is a rough reading of the observed baselines in Figure 1 (not a target). Neither point is named by the paper." * STAGE_NOTE,
        basis = :assumption, stage = :after_comparison, figures = ["fig4"]),
    Fix("mouse_dose_text_guess", :missing,
        "Mouse panel simulated with the middle dose group, 30 mg/kg (AUC 656.8 ug*h/mL).";
        evidence = "Eq. 2-3 text gives AUC for 10, 30 and 100 mg/kg (220.2, 656.8, 1140.8); the caption does not say which one Figure 4 left shows. The middle group is an arbitrary choice among the three." * STAGE_NOTE,
        basis = :assumption, stage = :after_comparison, figures = ["fig4"]),
    PATIENT_AUC_FIX, DAYS_PER_MONTH_FIX, CONTINUOUS_FIX,
]

const FIXES = [
    Fix("p0", :from_target,
        "Initial tumour size P0 is not stated for the typical simulation; set to 150 mm^3 (mice) and 66.8 mm (patients), read from Figure 4 at t = 0.";
        evidence = "Methods: P0 'set to the observed value' (individual baselines); no typical value given in text, Table 1, s001.txt/s004.txt (P0 = BASE/VB data column). Value read from the curve start in PSP4-4-660-g004.jpg (circular).",
        basis = :figure, stage = :after_comparison, figures = ["fig4"]),
    Fix("mouse_auc", :from_target,
        "Mouse panel simulated with the 100 mg/kg exposure (AUC 1140.8 ug*h/mL); the paper does not say which dose group the typical mouse simulation uses.";
        evidence = "Only gamma depends on AUC (beta_alpha fixed to 0). With P0 = 150: 100 mg/kg rms 0.001, 30 mg/kg rms 0.074, 10 mg/kg rms 0.27 (explore log). Chosen by matching the figure.",
        basis = :figure, stage = :after_comparison, figures = ["fig4"]),
    PATIENT_AUC_FIX, DAYS_PER_MONTH_FIX, CONTINUOUS_FIX,
]

function setting(ref_id::AbstractString, fixes)
    ids = Set(f.id for f in fixes)
    mouse = startswith(ref_id, "fig4_left")
    P0 = "p0" in ids ? (mouse ? FIXED_P0.mouse : FIXED_P0.patient) :
         "p0_text_guess" in ids ? (mouse ? TEXT_ONLY.mouse.P0 : TEXT_ONLY.patient.P0) :
         throw(NotRunnable("the typical initial tumour size P0 is not stated (Methods: 'set to the observed value')"))
    AUC = if mouse
        "mouse_auc" in ids ? FIXED_MOUSE_AUC : "mouse_dose_text_guess" in ids ? TEXT_ONLY.mouse.AUC :
            throw(NotRunnable("the dose group of the typical mouse simulation is not stated"))
    else
        "patient_auc_mean" in ids ? PATIENT_AUC : throw(NotRunnable("the exposure of the typical patient simulation is not stated"))
    end
    mouse || "days_per_month" in ids || throw(NotRunnable("the days-to-months conversion of Figure 4 right is not stated"))
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
