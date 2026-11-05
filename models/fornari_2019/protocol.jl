# Protocol for the fornari_2019 reproduction: dosing, simulation, frozen targets and tolerances, fixes, benchmark setups.
#
# TARGETS FROZEN 2026-11-04 09:49 IST (reference/*.csv written by reference/digitise_run.jl, tolerance below fixed, before any
# simulation of this model was run; see notes.md).
#
# Graded targets (red "mean model prediction" lines digitised from the figure images):
#   Figure 3 (rat, carboplatin 40 mg/kg iv on day 1 of eight 14-day cycles): MPP, CMP, MEP (bone marrow), neutrophils, platelets (blood)
#   Figure 4 (human, two 21-day carboplatin cycles, 42 days): neutrophils and platelets, panel (a) rat slope values, panel (b) Table 2 human slopes
# Not graded (listed for completeness): Fig 3 monocytes, reticulocytes, red blood cells and all observation symbols/percentile bands,
#   Fig 4 blue Schmitt-model lines, percentile bands and the bar chart (c), Supplement Figures S1, S2, S4-S7.

using SciMLBase: solve, successful_retcode
using OrdinaryDiffEq: Rodas5P

const TOLERANCE = Tolerance()   # default: RMS 0.015, P95 0.03 of the axis box; no per-curve overrides
const TOLERANCES = Dict{String, Tolerance}()

# ---------------------------------------------------------------------------------------------------------------------
# Fixes (each needs evidence; none re-estimates a parameter)
# ---------------------------------------------------------------------------------------------------------------------
# Provenance ledger (2026-11-05): every unstated or corrected choice is a Fix with basis and stage. The as-published run
# needs no fill: it stops at the unclosed equations (NotRunnable) before any dose or human value is needed.
const PUBLISHED_FILLS = Fix[]

const PRE = " Stage: listed in notes.md before any simulation was compared (equations, 'Missing / ambiguous inputs found', 'Human-parameter rule')."
const FIXES = [
    Fix("eq1_maturation", :typo,
        "Eq 1 maturation terms read -ktr_CMP fdbk_mat_MPP - ktr_MEP fdbk_mat_MPP. Used: -ktr_CMP fdbk_mat_CMP MPP - ktr_MEP fdbk_mat_MEP MPP.";
        evidence = "fdbk_mat_MPP is not defined anywhere (Eq 12 defines fdbk_mat_CMP, fdbk_mat_MEP); Eq 4 gain terms ktr_CMP fdbk_mat_CMP MPP and ktr_MEP fdbk_mat_MEP MPP are the same fluxes; supplement S4 kstem = (ktrCMP+ktrMEP-kprolMPP) MPP0 only closes the MPP steady state if the outflow is proportional to MPP; text on Eq 12 says the feedbacks are 'integrated in Eqs 1 and 4'." * PRE,
        basis = :paper, stage = :before_comparison),
    Fix("rbc_equation", :missing,
        "Eq 9 prints only Neut, Mono, Ret, Plt. Added dRBC/dt = kcirc_Ret Ret - kcirc_RBC RBC.";
        evidence = "Table 1 lists kcirc_RBC and the RBC baseline; Eqs 12-13 use RBC as a feedback signal; supplement S4 defines kcirc_Ret = kcirc_RBC RBC0/Ret0 (steady state of exactly this equation); S2 lists RBC among the states." * PRE,
        basis = :supplement, stage = :before_comparison),
    Fix("plt_transit_rate", :typo,
        "Eq 7 platelet T1 equation ends with -a_Ret T1_Plt. Used -a_Plt T1_Plt.";
        evidence = "Eq 8 defines a_Plt; the T2_Plt equation and supplement S3/S4 (ktrPlt = aPlt T1Plt/MEP0 ...) use a_Plt for the platelet chain; with a_Ret the S3 steady state of T1Plt is violated." * PRE,
        basis = :paper, stage = :before_comparison),
    Fix("typical_value_lines", :missing,
        "The red 'average model prediction' lines of Fig 3 and Fig 4 are simulated as a single typical-value run.";
        evidence = "The lines are averages of 1000 simulated subjects, but no inter-individual variability model is given, so only a typical-value simulation can be run. The Table S4 residual error shifts the mean by about 1.5-1.8 % (report.md, Issues 14)." * PRE,
        basis = :assumption, stage = :before_comparison),
    Fix("dose_time", :ambiguity,
        "Rat doses given at t = 24 h + 14 d x (0..7) (dose on 'day 1' of the time axis, which starts at day 0).";
        evidence = "Supplement Fig S3 plots the 30 mg/kg PK and Damage on an hour axis starting at 24 h with the dose at the left edge; Table S1 'Carbo dosing days 1, 15, 29...'; the text says eight cycles (Table S1 lists 1,15,...,85 for study E in the PDF)." * PRE,
        basis = :supplement, stage = :before_comparison, figures = ["fig3"]),
    Fix("carboplatin_mw", :missing,
        "Carboplatin molecular weight 371.25 g/mol converts mg doses to umol (rat: 2.694 umol per mg; human: 750 mg = 2020 umol).";
        evidence = "The paper gives no molecular weight; 371.25 g/mol is the standard value (PubChem CID 426756)." * PRE,
        basis = :standard, stage = :before_comparison),
    Fix("pk_per_kg", :units,
        "Table S3 rat PK printed as QC 0.06 '1/h', VCen 0.26 'l', VPer 0.30 'l', CL 0.42 '1/h' (code: liter/hour, liter). Read per kg (L/kg, L/h/kg) with the dose per kg: d mg/kg enters Cen as d x 1000 / 371.25 umol/kg (2.694 umol/kg per mg/kg; 107.7 at 40 mg/kg), Cp = Cen/VCen in umol/L.";
        evidence = "(1) Table S3 was estimated by fitting Eq S1 to the Fig S3 free-Pt data (Table S3 caption). Read per kg, the model reproduces those data (in-house 30 mg/kg, 16 points: geometric mean obs/model 0.93; Siddik 20 mg/kg, 10 points: 1.22); read as litres for a 0.320 kg rat it is 2.9x and 3.8x below them. These are PK data, not targets. " *
                   "(2) The cited PK source, Siddik et al. 1987 (Biochem Pharmacol 36:1925; abstract only, full text not accessed), reports carboplatin plasma clearance 10.3 mL/min/kg and inulin clearance 10.1 mL/min/kg in a study of mice and rats (the abstract does not say which species these values are from); read per kg, CL 0.42 L/h/kg is 7.0 mL/min/kg, the same order. " *
                   "As absolute litres for a 0.320 kg rat (S1) CL is 21.9 mL/min/kg, implausibly twice the GFR of a drug cleared mainly by filtration, and VCen 0.81 L/kg (Vss 1.75 L/kg) exceeds total body water (~0.67 L/kg) for a hydrophilic drug; per kg VCen 0.26 L/kg is extracellular water. " *
                   "(3) The Fig S3 data also give a model-free clearance of 0.47 L/h per kg (0.15 L/h for a 0.320 kg rat) against the printed 0.42. " *
                   "No value is fitted: Table S3 numbers are used unchanged; the 0.320 kg absolute reading (0.862 umol per mg/kg) is a sensitivity run. Nothing in the paper says 'per kg': the fitted values behave as per-kg values. " *
                   "Stage: adopted after an independent rebuild found that this reading matches Fig 3; the evidence above does not use Fig 3 or Fig 4 (selfcheck units resolution, 2026-11-05).",
        basis = :supplement, stage = :after_comparison, figures = ["fig3"]),
    Fix("human_schedule", :ambiguity,
        "Human simulation: two carboplatin doses of 750 mg, at day 0 and day 21 (Calvert dose AUC x (GFR + 25) = 5 x (125 + 25) mg, S11).";
        evidence = "Fig 4 caption (main text): 'Two cycles of 21 days of carboplatin were simulated. Doses were calculated with the Calvert formula ... targeting an AUC ...'. Supplement S11/S12 say 'a single dose'; the figure-specific caption is used because it describes the plotted simulation. " *
                   "Stage: added in the review round after the first comparison, on the caption text (not chosen by matching).",
        basis = :paper, stage = :after_comparison, figures = ["fig4"]),
    Fix("human_pk", :missing,
        "Human PK from Zandvliet 2008: CL = 0.76 x CrCl + 1.5 L/h, V1 15.5 L, Q 3.46 L/h, V2 9.86 L; each dose as a 30 min infusion (as in Zandvliet); ultrafiltrable Pt in umol/L taken as carboplatin umol/L.";
        evidence = "Supplement S11 states the human PK model is the 2-compartment model of Zandvliet et al. (Br J Clin Pharmacol 2008;66:485-497, PMC2561111); the parameter values and the infusion duration are in Zandvliet (Table 2, study design), not in the paper or supplement." * PRE,
        basis = :reference, stage = :before_comparison, figures = ["fig4"]),
    Fix("human_gfr_as_crcl", :missing,
        "The S11 GFR of 125 mL/min (7.5 L/h) is used as the creatinine clearance in the Zandvliet CL formula: CL = 0.76 x 7.5 + 1.5 = 7.2 L/h.";
        evidence = "S11 gives GFR = 125 mL/min for the Calvert dose; Zandvliet's CL is a function of Cockcroft-Gault creatinine clearance. The paper does not link the two." * PRE,
        basis = :assumption, stage = :before_comparison, figures = ["fig4"]),
    Fix("human_baselines", :ambiguity,
        "Table 1 prints human baselines for Neut, Mono, Ret, RBC as ranges (2-7, 0.2-10, 40-115, 4100-5900); the midpoints 4.5, 5.1, 77.5, 5000 are used.";
        evidence = "Table 1 human column. The paper names no point in the ranges, so the midpoint is an assumption. The midpoint rule was declared in notes.md before any simulation and is not read from Figure 4 (the figure shows ~4.7 neutrophils at t = 0). The value matters: Mono0 changes ktr_Mono and thus the CMP/neutrophil time scale." * PRE,
        basis = :assumption, stage = :before_comparison, figures = ["fig4"]),
]

# Sensitivity only (not graded): Table S3 read as absolute litres for a 0.320 kg rat (S1 minimum weight).
const PK_ABSOLUTE_SENSITIVITY = Fix("pk_absolute_320g", :units,
    "Sensitivity only: Table S3 read as absolute litres for a 0.320 kg rat, dose 0.320 kg x 1000/371.25 = 0.862 umol per mg/kg.";
    evidence = "S1: rats of minimum 320 g. Used in sensitivity.jl only.", basis = :assumption, stage = :after_comparison, figures = ["fig3"])

fixhas(fixes, id) = any(f -> f.id == id, fixes)

const DOSE_UMOL_PER_MGKG_S3 = 117.0 / 30               # diagnostic only (dose_scale_from_S3.jl, sensitivity.jl), never used in reproduce()
const DOSE_UMOL_PER_MGKG_PERKG = 1000 / 371.25         # fix pk_per_kg: umol/kg per mg/kg into VCen 0.26 L/kg (MW 371.25 g/mol)
const DOSE_UMOL_PER_MGKG_PHYS = 0.320 * 1000 / 371.25   # sensitivity only: Table S3 as absolute litres, 0.320 kg rat (S1 minimum weight)
const HUMAN_DOSE_UMOL = 750.0 / 371.25 * 1000          # Calvert: 5 x (125 + 25) = 750 mg
const HUMAN_INFUSION_H = 0.5

# ---------------------------------------------------------------------------------------------------------------------
# Problems
# ---------------------------------------------------------------------------------------------------------------------
function rate_callback(sys, param, times, values)
    ts = collect(Float64, times); vs = collect(Float64, values)
    set! = QSPAtlas.setu(sys, param)
    affect!(integrator) = begin
        i = findfirst(x -> isapprox(x, integrator.t; atol = 1e-9, rtol = 1e-12), ts)
        i === nothing || set!(integrator, vs[i])
    end
    return QSPAtlas.PresetTimeCallback(ts, affect!; save_positions = (true, true))
end

"""Rat problem: `ndoses` doses of `dose_mgkg` mg/kg every 14 days. Returns (prob, sys, callback, saveat, tend_hours)."""
function rat_problem(fixes; dose_mgkg = 40.0, ndoses = 8, tend_d = 120.0, dose_scale = nothing, rat_drug_params = rat_drug(), dense = 0.05)
    sys, u0 = haematopoiesis_system(rat_physiology(), rat_drug_params, rat_pk();
        structure = (fixhas(fixes, "eq1_maturation") && fixhas(fixes, "rbc_equation") && fixhas(fixes, "plt_transit_rate")) ? :corrected : :as_printed)
    scale = if dose_scale !== nothing
        dose_scale
    else
        fixhas(fixes, "carboplatin_mw") || throw(QSPAtlas.NotRunnable("the molecular weight for the mg/kg to umol conversion is not stated"))
        fixhas(fixes, "pk_per_kg") ? DOSE_UMOL_PER_MGKG_PERKG : fixhas(fixes, "pk_absolute_320g") ? DOSE_UMOL_PER_MGKG_PHYS :
            throw(QSPAtlas.NotRunnable("the dose-to-amount conversion (body weight or per-kg basis) is not stated"))
    end
    t0 = fixhas(fixes, "dose_time") ? 24.0 : 0.0
    times = repeated_doses(t0, 336.0, ndoses)
    cb = bolus_callback(sys, sys.Cen, times, dose_mgkg * scale)
    saveat = collect(0.0:(dense * 24):(tend_d * 24))
    prob = ODEProblem(sys, u0, (0.0, tend_d * 24); callback = cb)
    return (; prob, sys, callback = cb, saveat, tend = tend_d * 24)
end

"""Human problem (two carboplatin doses with fix human_schedule, else one). `slopes` is :rat (Figure 4a) or :table2 (Figure 4b) or :eq10 (sensitivity only)."""
function human_problem(fixes; slopes = :table2, tend_d = 42.4, phys = nothing, pk = nothing, infusion = HUMAN_INFUSION_H, dose_umol = HUMAN_DOSE_UMOL, dense = 0.05, dose_times_h = nothing)
    for id in ("eq1_maturation", "rbc_equation", "plt_transit_rate", "carboplatin_mw", "human_pk", "human_gfr_as_crcl", "human_baselines")
        fixhas(fixes, id) || throw(QSPAtlas.NotRunnable("human model needs fix $id (human PK and several human baselines are not given in the paper)"))
    end
    dose_times_h === nothing && (dose_times_h = fixhas(fixes, "human_schedule") ? [0.0, 504.0] : [0.0])
    drug = slopes === :rat ? rat_drug() : slopes === :table2 ? human_drug_table2() : human_drug_eq10()
    sys, u0 = haematopoiesis_system(phys === nothing ? human_physiology() : phys, drug, pk === nothing ? human_pk() : pk)
    rate = infusion > 0 ? dose_umol / infusion : 0.0
    if infusion > 0
        # first dose at t = 0 starts through the initial value of `rin`; every later start and every stop is a preset-time event
        ts = Float64[]; vs = Float64[]
        for td in dose_times_h
            td > 0 && (push!(ts, td); push!(vs, rate))
            push!(ts, td + infusion); push!(vs, 0.0)
        end
        p = sortperm(ts)
        cb = rate_callback(sys, sys.rin, ts[p], vs[p])
        u = vcat(u0, [sys.rin => rate])
    else
        cb = bolus_callback(sys, sys.Cen, dose_times_h, dose_umol); u = u0
    end
    saveat = collect(0.0:(dense * 24):(tend_d * 24))
    prob = ODEProblem(sys, u, (0.0, tend_d * 24); callback = cb)
    return (; prob, sys, callback = cb, saveat, tend = tend_d * 24)
end

const _CACHE = Dict{Any, Any}()
function solved(key, build)
    get!(_CACHE, key) do
        s = build()
        sol = solve(s.prob, Rodas5P(); abstol = 1e-9, reltol = 1e-8, saveat = s.saveat)
        (; s..., sol, ok = successful_retcode(sol))
    end
end

fixkey(fixes) = join(sort([f.id for f in fixes]), ",")

"""Simulated solution for a case: :rat (40 mg/kg x8), :human_a (rat slopes), :human_b (Table 2 slopes)."""
function case_solution(case::Symbol, fixes)
    key = (case, fixkey(fixes))
    case === :rat && return solved(key, () -> rat_problem(fixes))
    case === :human_a && return solved(key, () -> human_problem(fixes; slopes = :rat))
    case === :human_b && return solved(key, () -> human_problem(fixes; slopes = :table2))
    error("unknown case $case")
end

# reference id -> (case, variable)
const TARGET_MAP = Dict(
    "fig3_MPP" => (:rat, :MPP), "fig3_CMP" => (:rat, :CMP), "fig3_MEP" => (:rat, :MEP), "fig3_Neut" => (:rat, :Neut), "fig3_Plt" => (:rat, :Plt),
    "fig4a_neut" => (:human_a, :Neut), "fig4a_plt" => (:human_a, :Plt), "fig4b_neut" => (:human_b, :Neut), "fig4b_plt" => (:human_b, :Plt),
)

"""Time (days) and value of state `v` for a case."""
function series(case::Symbol, v::Symbol, fixes)
    s = case_solution(case, fixes)
    x = s.sol.t ./ 24
    # a failed solve gives NaN values (scored as a gap), never a partial trajectory
    s.ok || return x, fill(NaN, length(x))
    y = getproperty(s.sys, v)
    return x, [u for u in s.sol[y]]
end

simulate_curve(ref; fixes = Fix[]) = begin
    case, v = TARGET_MAP[ref.id]
    series(case, v, fixes)
end

# ---------------------------------------------------------------------------------------------------------------------
# Solver benchmark problems
# ---------------------------------------------------------------------------------------------------------------------
function benchmark_setups()
    r = rat_problem(FIXES)
    h = human_problem(FIXES; slopes = :table2)
    return [
        BenchmarkSetup(; name = "fornari_2019 rat, carboplatin 40 mg/kg x8 q14d (Fig 3)", prob = r.prob, saveat = r.saveat),
        BenchmarkSetup(; name = "fornari_2019 human, two 21-day cycles, Table 2 slopes (Fig 4b)", prob = h.prob, saveat = h.saveat),
    ]
end
