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

using SciMLBase: solve
using OrdinaryDiffEq: Rodas5P

const TOLERANCE = Tolerance()   # default: RMS 0.015, P95 0.03 of the axis box; no per-curve overrides
const TOLERANCES = Dict{String, Tolerance}()

# ---------------------------------------------------------------------------------------------------------------------
# Fixes (each needs evidence; none re-estimates a parameter)
# ---------------------------------------------------------------------------------------------------------------------
const FIXES = [
    Fix("eq1_maturation", :typo,
        "Eq 1 maturation terms read -ktr_CMP fdbk_mat_MPP - ktr_MEP fdbk_mat_MPP. Used: -ktr_CMP fdbk_mat_CMP MPP - ktr_MEP fdbk_mat_MEP MPP.";
        evidence = "fdbk_mat_MPP is not defined anywhere (Eq 12 defines fdbk_mat_CMP, fdbk_mat_MEP); Eq 4 gain terms ktr_CMP fdbk_mat_CMP MPP and ktr_MEP fdbk_mat_MEP MPP are the same fluxes; supplement S4 kstem = (ktrCMP+ktrMEP-kprolMPP) MPP0 only closes the MPP steady state if the outflow is proportional to MPP; text on Eq 12 says the feedbacks are 'integrated in Eqs 1 and 4'."),
    Fix("rbc_equation", :missing,
        "Eq 9 prints only Neut, Mono, Ret, Plt. Added dRBC/dt = kcirc_Ret Ret - kcirc_RBC RBC.";
        evidence = "Table 1 lists kcirc_RBC and the RBC baseline; Eqs 12-13 use RBC as a feedback signal; supplement S4 defines kcirc_Ret = kcirc_RBC RBC0/Ret0 (steady state of exactly this equation); S2 lists RBC among the states."),
    Fix("plt_transit_rate", :typo,
        "Eq 7 platelet T1 equation ends with -a_Ret T1_Plt. Used -a_Plt T1_Plt.";
        evidence = "Eq 8 defines a_Plt; the T2_Plt equation and supplement S3/S4 (ktrPlt = aPlt T1Plt/MEP0 ...) use a_Plt for the platelet chain; with a_Ret the S3 steady state of T1Plt is violated."),
    Fix("dose_time", :ambiguity,
        "Rat doses given at t = 24 h + 14 d x (0..7) (dose on 'day 1' of the time axis, which starts at day 0); human first dose at t = 0 (second at day 21, see human_schedule).";
        evidence = "Supplement Fig S3 plots the 30 mg/kg PK and Damage on an hour axis starting at 24 h with the dose at the left edge; Table S1 'Carbo dosing days 1, 15, 29...'; Fig 3 shows eight cycles (text: eight cycles) although Table S1 lists 1,15,...,85 for study E in the PDF. The Fig 4 axis starts at day 0 where the red line begins."),
    Fix("pk_per_kg", :units,
        "Table S3 rat PK printed as QC 0.06 '1/h', VCen 0.26 'l', VPer 0.30 'l', CL 0.42 '1/h' (code: liter/hour, liter). Read per kg (L/kg, L/h/kg) with the dose per kg: d mg/kg enters Cen as d x 1000 / 371.25 umol/kg (2.694 umol/kg per mg/kg; 107.7 at 40 mg/kg), Cp = Cen/VCen in umol/L.";
        evidence = "(1) The rat PK was fitted (S5, Table S3 caption) to in-house 30 mg/kg data and to Siddik et al. 1987 (Biochem Pharmacol 36:1925, supplement ref 1), which reports carboplatin plasma clearance 10.3 mL/min/kg, close to inulin 10.1 mL/min/kg (glomerular filtration only): per kg, CL 0.42 L/h/kg = 7.0 mL/min/kg is consistent; as absolute litres for a 0.320 kg rat (S1) it is 21.9 mL/min/kg, twice the GFR of a drug cleared by filtration alone, and VCen would be 0.81 L/kg (Vss 1.75 L/kg), more than total body water (~0.67 L/kg) for a hydrophilic drug; per kg VCen 0.26 L/kg is extracellular water. (2) The observed in-house 30 mg/kg free-Pt points of Supplement Fig S3 (not a target) at 1, 3, 6, 12, 24 h after dose are matched by the per-kg reading with mean log10 offset 0.00 (rms 0.15 decade); the absolute reading is 3.1x low (pk_units_check.jl). (3) All doses in the paper are per kg (S1: 'dose volume 5ml/kg', Table S1 mg/Kg). Carboplatin MW 371.25 g/mol (PubChem CID 426756). No value is fitted: Table S3 numbers are used unchanged; the 0.320 kg absolute reading (0.862 umol per mg/kg) is kept as a sensitivity run."),
    Fix("human_schedule", :ambiguity,
        "Human simulation: two carboplatin doses of 750 mg, at day 0 and day 21.";
        evidence = "Fig 4 caption (main text): 'Two cycles of 21 days of carboplatin were simulated. Doses were calculated with the Calvert formula ... targeting an AUC ...'. Supplement S11/S12 say 'a single dose'; the figure-specific caption is used because it describes the plotted simulation."),
    Fix("human_pk", :missing,
        "Human PK from Zandvliet 2008: CL = 0.76 x 7.5 + 1.5 = 7.2 L/h, V1 15.5 L, Q 3.46 L/h, V2 9.86 L; each 750 mg dose (= 2020 umol, MW 371.25) as a 30 min infusion.";
        evidence = "Supplement S11 states the human PK model is the 2-compartment model of Zandvliet et al. (Br J Clin Pharmacol 2008;66:485-497, PMC2561111) and the Calvert dose AUC (GFR+25) with AUC 5 mg min/mL, GFR 125 mL/min = 750 mg; the parameter values are in Table 2 of Zandvliet (not in the paper or supplement). Assumed: GFR = CG creatinine clearance, 30 min infusion as in Zandvliet, ultrafiltrable Pt in umol/L equals carboplatin umol."),
    Fix("human_baselines", :ambiguity,
        "Table 1 prints human baselines for Neut, Mono, Ret, RBC as ranges (2-7, 0.2-10, 40-115, 4100-5900); the midpoints 4.5, 5.1, 77.5, 5000 are used.";
        evidence = "Table 1 human column. The midpoint rule was declared in models/fornari_2019/notes.md before any simulation and is not read from Figure 4 (the figure shows ~4.7 neutrophils at t = 0). The value matters: Mono0 changes ktr_Mono and thus the CMP/neutrophil time scale."),
]

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
    scale = dose_scale !== nothing ? dose_scale : fixhas(fixes, "pk_per_kg") ? DOSE_UMOL_PER_MGKG_PERKG : DOSE_UMOL_PER_MGKG_PHYS
    t0 = fixhas(fixes, "dose_time") ? 24.0 : 0.0
    times = repeated_doses(t0, 336.0, ndoses)
    cb = bolus_callback(sys, sys.Cen, times, dose_mgkg * scale)
    saveat = collect(0.0:(dense * 24):(tend_d * 24))
    prob = ODEProblem(sys, u0, (0.0, tend_d * 24); callback = cb)
    return (; prob, sys, callback = cb, saveat, tend = tend_d * 24)
end

"""Human problem (two carboplatin doses with fix human_schedule, else one). `slopes` is :rat (Figure 4a) or :table2 (Figure 4b) or :eq10 (sensitivity only)."""
function human_problem(fixes; slopes = :table2, tend_d = 42.4, phys = nothing, pk = nothing, infusion = HUMAN_INFUSION_H, dose_umol = HUMAN_DOSE_UMOL, dense = 0.05, dose_times_h = nothing)
    for id in ("eq1_maturation", "rbc_equation", "plt_transit_rate", "human_pk", "human_baselines")
        fixhas(fixes, id) || error("human model needs fix $id (human PK and several human baselines are not given in the paper)")
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
        (; s..., sol)
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
        BenchmarkSetup(; name = "fornari_2019 rat, carboplatin 40 mg/kg x8 q14d (Fig 3)", prob = r.prob, callback = r.callback, saveat = r.saveat),
        BenchmarkSetup(; name = "fornari_2019 human, two 21-day cycles, Table 2 slopes (Fig 4b)", prob = h.prob, callback = h.callback, saveat = h.saveat),
    ]
end
