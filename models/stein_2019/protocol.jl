# Protocol, targets, tolerances, fixes for stein_2019 (Figure 3a, simulated percentile bands).
# TARGETS AND DESIGN FROZEN 2026-11-04 14:50 IST, before any simulation was compared with the
# bands (reference/*.csv, digitise.toml, the design constants and the tolerance below).
# No simulation had been run against the targets at that time (only the closed-form check of the
# model equations, which does not touch the bands). REVISED after review (see the design section and
# report.md): the as-published run now uses the scheduled visits; targets and tolerance unchanged except
# two misassigned points of fig3a_right_p90 removed by a documented rule. REVISED 2026-11-05 (after the
# first comparison): the "with fixes" run is the assumed VPC reconstruction informed by the Monolix documentation (RECON_FIXES below);
# the previous fixed design is kept as PREVIOUS_FIXES. Targets, reference CSVs and tolerance unchanged.
using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using Statistics: quantile, median

# ---- what is graded -----------------------------------------------------------------------
# Figure 3a ("Visual predictive check") has two panels (0-380 d and a zoom on 0-30 d). Blue lines are
# the 10/50/90th percentiles of the DATA (observed; not graded). The shaded bands are the
# authors' simulation: the prediction interval of the simulated 10th and 90th percentile (blue)
# and of the simulated 50th percentile (pink). Graded: where the simulated percentile of OUR
# simulation lies relative to those bands, one `:band` reference per band and panel.
# Default tolerance (RMS 0.015, P95 0.03 of the axis box): it is the allowance for digitising the
# band edges. With the band criterion an error is 0 inside the band, so the tolerance is not a
# fit tolerance on the percentile itself. Not changed after comparing.
const TOLERANCE = Tolerance()
const TOLERANCES = Dict{String, Tolerance}()

# ---- simulation design -----------------------------------------------------------------------
# REVISED after review (2026-11-04): the frozen design of 14:50 used a dense time grid
# and left-panel nodes read from the target for the graded run. The first (14:50) scores are in
# report.md. The graded runs are now:
#   "as published" = the most literal reading of Methods "Sample analysis": 90 patients (Methods),
#       each observed at the SCHEDULED visit days 1, 2, 7, 11, 14, 21, 28 and months 3 and 6
#       (91.3125 and 182.625 d; month = 365.25/12 d; later visits are beyond the figure), the
#       simulated percentile at those times joined by straight segments in (day, log10 copies).
#   "with fixes"   = until 2026-11-05 the same simulation read on a denser design: a dense grid (0.25 d on
#       (0, 30]) for the right panel and the left-panel bin nodes of the band (28, 68, 190 d) joined by
#       chords (PREVIOUS_FIXES). Since 2026-11-05 the assumed Monolix-informed VPC reconstruction with its own
#       simulation (see "VPC design reconstruction" below and RECON_FIXES).
# Choices the paper does not state and that BOTH runs need (they are in FIXES so they are visible,
# but the as-published run cannot be defined without them): every patient observed at every time
# (no missing visits, no unscheduled samples), residual error included, simulated BLQ values kept,
# quantile type, base model without comedication.
const N_PATIENTS = 90          # Methods: qPCR in 90 patients (ELIANA 61 + ENSIGN 29)
const REPLICATES = 800         # 400 gave a max MC SE of 0.0034 box, above the 0.003 rule set beforehand
const CHUNK = 50
const SEED = 20260411
const PI_LEVEL = 0.90          # level of the shaded interval in the regenerated figure only (Monolix default; the paper does not state it)
const QS = (0.1, 0.5, 0.9)
const DAYS_PER_MONTH = 365.25 / 12
const SCHEDULED = [1.0, 2.0, 7.0, 11.0, 14.0, 21.0, 28.0, 3 * DAYS_PER_MONTH, 6 * DAYS_PER_MONTH]
const TIMES_RIGHT = collect(0.25:0.25:30.0)
const TIMES_LEFT = collect(31.0:1.0:190.0)
const TIMES = sort(unique(vcat(TIMES_RIGHT, TIMES_LEFT, SCHEDULED)))
const NODES_LEFT = [28.0, 68.0, 190.0]   # bin nodes of the left panel, read from the band kinks (digitise.toml): :from_target
const TSPAN = (0.0, 190.0)

# ---- fixes (provenance ledger 2026-11-05) ---------------------------------------------------
# Every design choice that is not stated by the paper, with basis and stage. The shared fills are used by
# both runs (the as-published run is not definable without them); `scheduled_visits` and `days_per_month`
# define the as-published reading, `dense_grid` and `left_bin_nodes` the "with fixes" reading.
const FROZEN = " Stage: part of the design frozen at 14:50 before any comparison (notes.md)."
const SHARED_FILLS = [
    Fix("all_patients_all_times", :missing,
        "Every virtual patient is observed at every time of the design (no missing visits, no unscheduled samples).";
        evidence = "Methods 'Sample analysis' lists scheduled visits (days 1, 2, 7, 11, 14, 21, 28, then 3-monthly) and Results 836 measurements in 90 patients (9.3 per patient, 43 BLQ); who was sampled when is not given. Used by both runs." * FROZEN,
        basis = :assumption, stage = :before_comparison),
    Fix("residual_included", :ambiguity,
        "Percentiles are taken over observations WITH residual error (log-normal; Table 1 a = 0.56 read as the sd of the log-scale residual).";
        evidence = "Fig 3a plots observations against simulated percentiles of observations (a VPC); Table 1 residual error a = 0.56, Methods 'constant error model in log-space'. The caption does not say whether ipred or dv percentiles are shown. Used by both runs." * FROZEN,
        basis = :interpretation, stage = :before_comparison),
    Fix("blq_kept", :ambiguity,
        "Simulated values below the LLOQ (50 copies/ug) are kept, not replaced by the LLOQ.";
        evidence = "Caption: red asterisks are 'simulated data points for the measurements that were BLQ'. This shows simulated values are plotted; it does not say how they enter the percentiles. Censoring to the LLOQ is a sensitivity run. Used by both runs." * FROZEN,
        basis = :interpretation, stage = :before_comparison),
    Fix("quantile_type7", :missing,
        "Percentiles are Hyndman-Fan type 7 (linear interpolation of order statistics) over the 90 patients, median over replicates.";
        evidence = "Quantile convention of the VPC is not stated; type 7 is the R/Monolix-style default and is the QSPAtlas `percentiles` default. Used by both runs." * FROZEN,
        basis = :assumption, stage = :before_comparison),
    Fix("base_model_no_comedication", :ambiguity,
        "Base model (Eq. 1) without covariates and without the tocilizumab/corticosteroid effect (F = 1).";
        evidence = "Results 'Analytical model': 'The base model characterized the data well, as shown by the visual predictive check and individual fits in Figure 3'; Results 'Covariate analysis': the final model is the base model. Methods 'Covariate model' says the full covariate model is the final model, and Table 1 lists Ftoci = 1.2, so the VPC may have used comedication regressors; their per-patient values are not published. Comedication is a sensitivity run with assumed timings. Used by both runs." * FROZEN,
        basis = :interpretation, stage = :before_comparison),
    Fix("omega_as_sd", :ambiguity,
        "Table 1 random-effect numbers are read as standard deviations of eta (omega), not variances.";
        evidence = "Table 1 gives one unlabelled number per parameter; supplement line 77: 'random effect eta_i with variance omega^2'; Monolix reports omega as an SD. The variance reading is a sensitivity run. Used by both runs." * FROZEN,
        basis = :supplement, stage = :before_comparison),
    Fix("omega_diagonal", :missing,
        "Random effects are uncorrelated (diagonal Omega).";
        evidence = "No correlations between etas are reported; diagonal is the Monolix default. Used by both runs." * FROZEN,
        basis = :assumption, stage = :before_comparison),
]

const PUBLISHED_FILLS = [
    SHARED_FILLS...,
    Fix("scheduled_visits", :ambiguity,
        "As published: each patient observed at the scheduled visit days 1, 2, 7, 11, 14, 21, 28 and months 3 and 6 (day N read as t = N d); the simulated percentile at those times joined by straight segments in (day, log10 copies).";
        evidence = "Methods 'Sample analysis' lists the scheduled visits; this is the most literal reading of the design. Stage: adopted in the review round after the first comparison, as the literal reading (not chosen by matching; it scores slightly worse than the dense design).",
        basis = :interpretation, stage = :after_comparison),
    Fix("days_per_month", :missing,
        "Months 3 and 6 placed at 3 and 6 x 365.25/12 days (91.3 and 182.6 d).";
        evidence = "The paper gives the late visits in months; 365.25/12 is the mean calendar month. Stage: with scheduled_visits, after the first comparison.",
        basis = :standard, stage = :after_comparison),
]

# Previous "with fixes" design (2026-11-04 to 2026-11-05): dense right-panel grid and left-panel nodes read
# from the target. Kept for the record and the sensitivity runs; replaced as the graded fixed run by the
# assumed Monolix-informed reconstruction below (RECON_FIXES).
const PREVIOUS_FIXES = [
    SHARED_FILLS...,
    Fix("dense_grid", :missing,
        "Right panel read on a dense grid (0.25 d) instead of only at the scheduled visits.";
        evidence = "Monolix VPC bins come from the data's actual (unscheduled, irregular) sample times, which are not published; the band between visits is drawn from those bins. A dense grid removes the coarse chord between scheduled visits, which is not how the band was built. Only in the 'with fixes' run." * FROZEN,
        basis = :assumption, stage = :before_comparison),
    Fix("left_bin_nodes", :from_target,
        "Left panel simulated percentile taken at three bin nodes (28, 68, 190 d) and joined by chords.";
        evidence = "Nodes are read from the kinks of the authors' band edges in Figure 3a (digitise.toml); the bins themselves are not stated, so this input comes from the target figure (circular, the dagger case). Only in the 'with fixes' run." * FROZEN,
        basis = :figure, stage = :before_comparison),
]

# ---- VPC design reconstruction (2026-11-05) ---------------------------------------------------
# The graded "with fixes" run is an ASSUMED reconstruction of the VPC, informed by the Monolix documentation (not
# a verified copy of what the authors ran), from sources that are not
# the target figure: the paper's text, the supplement's mlxtran model, the Monolix documentation and
# Lavielle & Bleakley 2011 (the binning algorithm of Monolix 4). Chosen AFTER the first comparison (2026-11-04)
# and after the comparisons of the previous design, and the comedication regressors had been run as a sensitivity
# on 2026-11-04 and shown to improve agreement with the target, so the reconstruction is not blind to the target
# and its 4 of 6 bands are conditional on the assumptions. The paper does not state the Monolix version (the
# 2016R1 to 2019R1 documentation could not be retrieved): documented for 2024R1 are R5, 500 simulations, 90 % PIs,
# simulated BLQ in the percentiles; for Monolix 4 the least-squares binning; equal width is the 2024R1 default.
# So :reference only where a cited document states the setting; binning method, percentile convention, node
# position and regressor spreads are :assumption.
#   Monolix VPC (https://monolixsuite.slp-software.com/monolix/2024R1/visual-predictive-check): 10/50/90th
#     percentiles, 90 % prediction intervals, 500 simulations of the design of the data set, percentiles
#     "piecewise linear method, also called R5" (Hyndman-Fan type 5).
#   Censored data in the VPC (.../2024R1/correcting-the-vpc-for-missing-observations): with "use censored
#     data" on and BLQ "simulated", "the VPC is simulated for all times including censored times, and all
#     simulated values are used to calculate the simulated percentiles"; observed BLQ data are replaced by
#     draws from their conditional distribution (the red asterisks of the caption).
#   Binning: Lavielle M, Bleakley K. J Pharmacokinet Pharmacodyn 2011;38:861-871 (HAL hal-00654842),
#     implemented as the automatic binning of Monolix 4: bins minimise J = sum_k n_k (var_k)^beta (dynamic
#     programming, beta = 0.2), the number of bins minimises log(J) + lambda beta K (lambda = 0.3). The
#     2024R1 documentation lists least-squares binning as an option with equal width as the default and an
#     automatic bin count whose rule is not documented, so least-squares is the reproducible Monolix choice.
const MLX_PI_LEVEL = 0.90
const MLX_SIMULATIONS = 500           # Monolix default; the compared median uses RECON_REPLICATES (MC SE rule)
const MLX_BETA = 0.2
const MLX_LAMBDA = 0.3
const RECON_REPLICATES = 500          # MC SE rule (max < 0.003 box at every scored point) checked in reproduce()
const RECON_TSPAN = (0.0, 380.0)      # Results: "the longest follow-up for this analysis was only 1 year"
# Methods "Sample analysis": days 1, 2, 7, 11, 14, 21, 28, then 3-monthly from month 3 to 12 (follow-up 1 year)
const VISITS_12M = vcat([1.0, 2.0, 7.0, 11.0, 14.0, 21.0, 28.0], [3.0, 6.0, 9.0, 12.0] .* DAYS_PER_MONTH)

const RECON_STAGE = " Stage: VPC design reconstruction of 2026-11-05, after the first comparison (2026-11-04) and the review round; an assumed reconstruction informed by the sources named, after the 2026-11-04 comedication sensitivity had shown the regressors improve agreement with the bands; the 4 of 6 is conditional on the settings marked :assumption."
const RECON_FIXES = [
    filter(f -> f.id in ("all_patients_all_times", "residual_included", "omega_as_sd", "omega_diagonal"), SHARED_FILLS)...,
    Fix("mlx_base_model_regressors", :ambiguity,
        "Base model as in the supplement's mlxtran: Eq. 1 with the comedication effect of Eq. 2 entered per patient through the time-of-first-dose regressors (Ftoci 1.2, Fster 1.0), no covariate effects on Cmax.";
        evidence = "Results: 'The base model characterized the data well, as shown by the visual predictive check ... in Figure 3' and 'the final cellular kinetic model was the same as the base model'. The mlxtran printed in Supplementary Material S1 (lines 100-152, 'the mlxtran code used for fitting the model') has parameter = {foldx, Tmax, Cmax, Ftoci, Fster, alpha, FB, beta} and regressor = {TOCI1T, STERSTT}; Table 1 lists Ftoci and Fster as fixed effects. The covariate model of the paper is the set of Cmax covariates (Table 1 'log Cmax covariate effect'), not used here. A Monolix VPC simulates each subject with its own regressors. Replaces base_model_no_comedication of the previous design." * RECON_STAGE,
        basis = :supplement, stage = :after_comparison),
    Fix("mlx_regressor_timing", :missing,
        "Per-patient comedication regressors drawn from Table 2: tocilizumab with probability 0.36, first dose 5.7 d times exp(0.7 z) clamped to 1-27 d; corticosteroids with probability 0.26, 7.5 d times exp(1.2 z) clamped to 0.11-170 d (z standard normal); independent; none = no effect.";
        evidence = "The patients' own regressors are in Data set S1 (supplement, PSP4-8-285-s002.csv: 32 of 90 with tocilizumab, median 5.75 d; 23 with steroids, median 7.45 d), which is not redistributed with this repository, so the graded run cannot read them. Table 2 gives the proportions, medians and ranges; the log-normal shape, its spread and the independence are assumptions, the same sampler as the comedication sensitivity pre-registered on 2026-11-04, which had already shown that these regressors improve agreement with the bands (4 of 6), so this choice is not blind to the target; the spreads have no independent evidence. reproduce(s1 = true) runs the design of Data set S1 itself when the file is present locally." * RECON_STAGE,
        basis = :assumption, stage = :after_comparison),
    Fix("mlx_nominal_visits", :missing,
        "Each of the 90 patients observed at the Methods visits up to one year: days 1, 2, 7, 11, 14, 21, 28 and months 3, 6, 9, 12 (day N read as t = N d, month = 365.25/12 d).";
        evidence = "Methods 'Sample analysis' (visit schedule) and Results ('the longest follow-up for this analysis was only 1 year'). A Monolix VPC simulates the design of the data set (the actual sample times), which is in Data set S1 (not redistributed); the nominal visits are the published version of that design. Missing visits are not modelled (all_patients_all_times)." * RECON_STAGE,
        basis = :interpretation, stage = :after_comparison),
    Fix("days_per_month", :missing,
        "Months placed at multiples of 365.25/12 days.";
        evidence = "The paper gives the late visits in months; 365.25/12 is the mean calendar month." * RECON_STAGE,
        basis = :standard, stage = :after_comparison),
    Fix("mlx_ls_binning", :missing,
        "Bins of the pooled observation times by the Monolix least-squares criterion (Lavielle & Bleakley 2011: beta 0.2, number of bins by log(J) + 0.3 beta K). On the nominal visits every visit is its own bin.";
        evidence = "Lavielle & Bleakley, J Pharmacokinet Pharmacodyn 2011;38:861-871: the automatic binning implemented in Monolix 4 with these defaults. The paper does not state its Monolix version or bins; the Monolix 2024R1 documentation lists least-squares as an option with equal width as the default and an undocumented automatic bin count. With tied times a bin made of one visit has zero within-bin variance, so the criterion puts each visit in its own bin. Least-squares binning is therefore documented for Monolix 4 (and as an option in 2024R1, whose default is equal width), not for the version the paper used: an assumption." * RECON_STAGE,
        basis = :assumption, stage = :after_comparison),
    Fix("mlx_bin_node_mean_time", :missing,
        "Each bin's percentile is placed at the mean observation time of the bin and the nodes are joined by straight segments in (day, log10 copies).";
        evidence = "Lavielle & Bleakley 2011 take the bin mean time as the time that represents the bin; where Monolix draws the node is not documented. On the nominal visits the node is the visit time." * RECON_STAGE,
        basis = :assumption, stage = :after_comparison),
    Fix("mlx_percentile_r5", :missing,
        "Percentiles of each bin by Hyndman-Fan type 5 (piecewise linear, R5), not type 7.";
        evidence = "Monolix 2024R1 VPC documentation: 'The calculation of percentiles follows the piecewise linear method, also called R5 method.' Whether the version used for the paper did the same is not documented (the 2016R1 to 2019R1 pages could not be retrieved), so the convention is an assumption for the paper's era. QSPAtlas `percentiles` has type 7 only, so type 5 is computed in this model's protocol.jl." * RECON_STAGE,
        basis = :assumption, stage = :after_comparison),
    Fix("mlx_blq_simulated", :ambiguity,
        "All simulated values, including those below the LLOQ of 50 copies/ug, enter the simulated percentiles (no censoring of the simulated data).";
        evidence = "Methods: BLQ data were fixed-point censored (CENSORING column); caption: red asterisks are simulated values for the BLQ measurements (the Monolix 'BLQ simulated' display). Monolix documentation, 'Correcting the VPC for missing observations' (2024R1): with censored data used and simulated, 'all simulated values are used to calculate the simulated percentiles'. Documented for 2024R1; the paper's caption (simulated BLQ points shown) is consistent with it. Replaces blq_kept of the previous design (same numbers, now with a source)." * RECON_STAGE,
        basis = :reference, stage = :after_comparison),
    Fix("mlx_pi_90_500", :missing,
        "Prediction interval 90 % over the simulated replicates (Monolix default; 500 simulations). The compared quantity is the median over replicates of each simulated percentile, with as many replicates as the Monte Carlo rule needs.";
        evidence = "Monolix 2024R1 VPC documentation: 'By default 500 simulations are performed', 'By default, prediction intervals are computed with a level of 90%'; Lavielle & Bleakley 2011 (Monolix 4) also use 90 %. The interval level only sets the displayed band in figures/fig3.png; the score compares the median of each simulated percentile with the authors' band." * RECON_STAGE,
        basis = :reference, stage = :after_comparison),
]
const FIXES = RECON_FIXES

"Which reading of the band curves a fix list asks for (NotRunnable if the design is not given)."
function band_mode(fixes)
    ids = Set(f.id for f in fixes)
    recon = ("mlx_base_model_regressors", "mlx_regressor_timing", "mlx_nominal_visits", "days_per_month", "mlx_ls_binning",
        "mlx_bin_node_mean_time", "mlx_percentile_r5", "mlx_blq_simulated", "mlx_pi_90_500", "all_patients_all_times", "residual_included", "omega_as_sd", "omega_diagonal")
    all(in(ids), recon) && return :monolix
    all(f -> f.id in ids, SHARED_FILLS) || throw(NotRunnable("the VPC design (sampling, residual, BLQ, quantiles, model, random effects) is not stated"))
    ("dense_grid" in ids && "left_bin_nodes" in ids) && return :fixed
    ("scheduled_visits" in ids && "days_per_month" in ids) && return :scheduled
    throw(NotRunnable("the VPC times and bins are not stated"))
end

# ---- population model ---------------------------------------------------------------------
const ETAS = [:foldx, :Tmax, :Cmax, :α, :FB, :β]

"""
    population_model(; comedication = false, omega = :sd, solver = Vern9(), tol = 1e-10)

PopulationModel for the CAR-T model. `omega = :sd` reads Table 1's random-effect numbers as
standard deviations (default), `:variance` as variances (sensitivity only). `comedication = true`
turns on the tocilizumab / corticosteroid effect of Eq. 2 with per-subject covariates `toci`,
`ster` (day of first dose, `Inf` for none), as in the supplement's Monolix code.
"""
function population_model(; comedication::Bool = false, omega::Symbol = :sd, solver = Vern9(), tol = 1.0e-10,
        theta = TYPICAL, omegas = OMEGA, sigma = RESIDUAL_SD_LOG, tspan = TSPAN)
    sds = [omega === :sd ? omegas[n] : sqrt(omegas[n]) for n in ETAS]
    Ω = [sds[i]^2 * (i == j) for i in eachindex(ETAS), j in eachindex(ETAS)]
    function params(θ, η, c)
        foldx_i = θ.foldx * exp(η.foldx)
        Tmax_i = θ.Tmax * exp(η.Tmax)
        base = [foldx => foldx_i, Tmax => Tmax_i, α => θ.α * exp(η.α),
            FB => θ.FB * exp(η.FB), β => θ.β * exp(η.β)]
        if comedication
            # supplement lines 108-123: the first comedication is (t1, F1), the second (t2, F2);
            # no dose is 99999 in the data set, here Inf
            if c.toci <= c.ster
                t1, f1, t2, f2 = min(c.toci, Tmax_i), θ.Ftoci, min(c.ster, Tmax_i), θ.Fster
            else
                t1, f1, t2, f2 = min(c.ster, Tmax_i), θ.Fster, min(c.toci, Tmax_i), θ.Ftoci
            end
        else
            t1, f1, t2, f2 = Tmax_i, 1.0, Tmax_i, 1.0
        end
        return vcat(base, [T1 => t1, T2 => t2, F1 => f1, F2 => f2])
    end
    init = (θ, η, c) -> [E => θ.Cmax * exp(η.Cmax) / (θ.foldx * exp(η.foldx)), M => 0.0]   # R0 = Cmax / foldx
    return PopulationModel(CART; routes = Route[], θ = theta, ηnames = ETAS, Ω, params, init,
        observed = (C = CART.C,), error = (C = ExponentialError(sigma),), tspan,
        solver, solve_kwargs = (abstol = tol, reltol = tol))
end

# The CAR-T "dose" is not a model input: the paper found no dose-exposure relationship and the
# infused cells enter only through the initial value R0 = Cmax / foldx. Empty regimen.
function cohort(n = N_PATIENTS; times = TIMES, comedication = nothing)
    cov = comedication === nothing ? (rng -> (toci = Inf, ster = Inf)) : comedication
    return virtual_population(n; regimen = Regimen(Dose[]), obstimes = times, covariate_sampler = cov, seed = SEED)
end

# ---- closed form (Eq. 1 / Eq. 2 and the supplement's Monolix code), for the model check ---------
function closed_form(t_, p)   # p: (foldx, Tmax, Cmax, α, FB, β, F1, F2, t1, t2)
    ρ_ = log(p.foldx) / p.Tmax
    R0 = p.Cmax / p.foldx
    R1 = R0 * exp(ρ_ * p.t1)
    R2 = R1 * exp(p.F1 * ρ_ * (p.t2 - p.t1))
    Cadj = R2 * exp(p.F1 * p.F2 * ρ_ * (p.Tmax - p.t2))
    t_ < p.t1 && return R0 * exp(ρ_ * t_)
    t_ < p.t2 && return R1 * exp(p.F1 * ρ_ * (t_ - p.t1))
    t_ < p.Tmax && return R2 * exp(p.F1 * p.F2 * ρ_ * (t_ - p.t2))
    return (1 - p.FB) * Cadj * exp(-p.α * (t_ - p.Tmax)) + p.FB * Cadj * exp(-p.β * (t_ - p.Tmax))
end

# ---- replicate percentiles in chunks --------------------------------------------------------
"""
    replicate_percentiles(model, subjects; replicates, chunk, seed, q, blq) -> (times, Q)

`Q[i, j, k]` is the q_j-quantile (type 7) of the dv of the cohort at time i in replicate k.
Cohort replicates are simulated in blocks of `chunk` (block b uses seed `seed + b`) to keep memory
bounded; the result does not depend on the thread count.
"""
function replicate_percentiles(model, subjects; replicates = REPLICATES, chunk = CHUNK, seed = SEED, q = QS, blq = nothing)
    times = sort(unique(subjects[1].obstimes))
    Q = fill(NaN, length(times), length(q), replicates)
    k0 = 0
    b = 0
    while k0 < replicates
        m = min(chunk, replicates - k0)
        res = simulate(model, subjects; replicates = m, seed = seed + b)
        nbad = count(s -> s.retcode !== :Success, res.subjects)
        nbad == 0 || error("stein_2019: $nbad of $(length(res.subjects)) simulated subjects did not finish with retcode Success in replicate block $b")
        rows = percentiles(res, :C; q, bins = :obstimes, value = :dv, by_replicate = true, blq)
        for r in rows
            i = searchsortedfirst(times, r.t_mid)
            Q[i, findfirst(==(r.q), collect(q)), k0 + r.replicate] = r.value
        end
        k0 += m
        b += 1
    end
    return times, Q
end

"""
    summarise(times, Q; pi_level) -> NamedTuple

Per time and quantile: median over replicates (the compared quantity), the central `pi_level`
interval across replicates, and the Monte Carlo standard error of the median (bootstrap over replicates, as in `vpc`).
"""
function summarise(times, Q; pi_level = PI_LEVEL, seed = SEED)
    nt, nq, K = size(Q)
    med = [median(Q[i, j, :]) for i in 1:nt, j in 1:nq]
    lo = [quantile(Q[i, j, :], (1 - pi_level) / 2) for i in 1:nt, j in 1:nq]
    hi = [quantile(Q[i, j, :], (1 + pi_level) / 2) for i in 1:nt, j in 1:nq]
    # same bootstrap estimator and seeding as QSPAtlas `vpc` (Part 4), so the numbers agree with `vpc`
    se = [QSPAtlas.mc_se_median(Q[i, j, :], QSPAtlas.stream_rng(seed, i, string(QS[j]), :mcse)) for i in 1:nt, j in 1:nq]
    return (; times, median = med, lo, hi, mc_se = se, replicates = K)
end

# ---- assumed Monolix-informed VPC (reconstruction) ---------------------------------------------------
# Comedication regressors drawn from Table 2 (see mlx_regressor_timing; the same sampler as the comedication
# sensitivity run pre-registered on 2026-11-04): per patient the day of the first tocilizumab and corticosteroid
# dose, Inf for none (the mlxtran uses 99999, which is clipped at Tmax all the same).
function comedication_sampler(rng)
    draw(median_, σ, lo, hi) = clamp(median_ * exp(σ * randn(rng)), lo, hi)
    toci = rand(rng) < 0.36 ? draw(5.7, 0.7, 1.0, 27.0) : Inf
    ster = rand(rng) < 0.26 ? draw(7.5, 1.2, 0.11, 170.0) : Inf
    return (toci = toci, ster = ster)
end

"""
    ls_binning(times; β = MLX_BETA, λ = MLX_LAMBDA, kmax = 100) -> Vector{NamedTuple}

Least-squares binning of Lavielle & Bleakley (2011), the automatic binning of Monolix 4. The pooled
observation times are cut into contiguous bins of distinct time values that minimise
J = sum_k n_k (var_k)^β (var_k the within-bin variance of the times, n_k the number of observations), by
dynamic programming, and the number of bins K minimises log(J) + λ β K. A tied time is never split. If J
reaches 0 (every distinct time its own bin, e.g. nominal visits) that K is taken (log 0 = -Inf).
Returns `(lo, hi, mid, n)` per bin; `mid` is the mean observation time of the bin.
"""
function ls_binning(times::AbstractVector{<:Real}; β = MLX_BETA, λ = MLX_LAMBDA, kmax::Integer = 100)
    t = sort(collect(Float64, times))
    u = unique(t)
    w = [searchsortedlast(t, x) - searchsortedfirst(t, x) + 1 for x in u]
    M = length(u)
    cw = cumsum(vcat(0, w)); c1 = cumsum(vcat(0.0, w .* u)); c2 = cumsum(vcat(0.0, w .* u .^ 2))
    function cost(i, j)                       # bin made of distinct times i..j
        n = cw[j + 1] - cw[i]
        m1 = (c1[j + 1] - c1[i]) / n
        v = max((c2[j + 1] - c2[i]) / n - m1^2, 0.0)
        v < 1.0e-12 * max(1.0, m1^2) && (v = 0.0)  # rounding of a single tied time
        return n * v^β
    end
    K = min(kmax, M)
    F = fill(Inf, K, M); B = zeros(Int, K, M)
    for j in 1:M
        F[1, j] = cost(1, j)
    end
    for k in 2:K, j in k:M, i in k:j
        c = F[k - 1, i - 1] + cost(i, j)
        if c < F[k, j]
            F[k, j] = c; B[k, j] = i
        end
    end
    U = [log(F[k, M]) + λ * β * k for k in 1:K]
    kbest = argmin(U)
    starts = Int[]; j = M
    for k in kbest:-1:2
        i = B[k, j]; pushfirst!(starts, i); j = i - 1
    end
    pushfirst!(starts, 1)
    ends = vcat(starts[2:end] .- 1, M)
    return [(lo = u[a], hi = u[b], mid = (c1[b + 1] - c1[a]) / (cw[b + 1] - cw[a]), n = cw[b + 1] - cw[a]) for (a, b) in zip(starts, ends)]
end

"Hyndman-Fan type 5 quantile (Monolix 'R5', piecewise linear) of an already sorted vector."
function quantile_type5(x::AbstractVector{Float64}, p::Real)
    n = length(x)
    n == 0 && return NaN
    h = n * p + 0.5
    h <= 1 && return x[1]
    h >= n && return x[n]
    f = floor(Int, h)
    return x[f] + (h - f) * (x[f + 1] - x[f])
end

"""
    monolix_percentiles(model, subjects, bins; replicates, chunk, seed) -> Q

`Q[k, j, r]` is the type-5 `QS[j]` percentile of ALL simulated observations (residual error included,
values below the LLOQ kept, Monolix 'BLQ simulated') that fall in bin `k` in replicate `r`. Simulated in
blocks of `chunk` replicates with seed `seed + block`, as `replicate_percentiles`.
"""
function monolix_percentiles(model, subjects, bins; replicates = RECON_REPLICATES, chunk = CHUNK, seed = SEED)
    nb = length(bins)
    Q = fill(NaN, nb, length(QS), replicates)
    binof(t_) = findfirst(b -> b.lo <= t_ <= b.hi, bins)
    k0 = 0; blk = 0
    while k0 < replicates
        m = min(chunk, replicates - k0)
        res = simulate(model, subjects; replicates = m, seed = seed + blk)
        nbad = count(s -> s.retcode !== :Success, res.subjects)
        nbad == 0 || error("stein_2019: $nbad simulated subjects did not finish with retcode Success in replicate block $blk")
        cells = [Float64[] for _ in 1:nb, _ in 1:m]
        for o in res.obs
            k = binof(o.time)
            k === nothing && error("stein_2019: observation at t = $(o.time) is in no bin")
            push!(cells[k, o.replicate], o.dv)
        end
        for k in 1:nb, r in 1:m
            v = sort!(cells[k, r])
            for j in eachindex(QS)
                Q[k, j, k0 + r] = quantile_type5(v, QS[j])
            end
        end
        k0 += m; blk += 1
    end
    return Q
end

# The reconstructed (graded) design: 90 patients at the nominal visits up to one year, regressors from Table 2.
recon_cohort(n = N_PATIENTS) = cohort(n; times = VISITS_12M, comedication = comedication_sampler)
recon_model(; kw...) = population_model(; comedication = true, tspan = RECON_TSPAN, kw...)

"Summary of a Monolix-style VPC: bins from the pooled design times, type-5 percentiles, 90 % interval."
function monolix_summary(model, subjects; replicates = RECON_REPLICATES, seed = SEED, bins = ls_binning(reduce(vcat, [s.obstimes for s in subjects])))
    Q = monolix_percentiles(model, subjects, bins; replicates, seed)
    S = summarise([b.mid for b in bins], Q; pi_level = MLX_PI_LEVEL, seed)
    return merge(S, (; bins))
end

# curve of one band target from a summary. mode :scheduled = scheduled visits joined by chords (as
# published); :monolix = the bin nodes of a Monolix-style summary (with fixes, reconstruction of 2026-11-05);
# :fixed = right panel dense, left panel bin nodes joined by chords (previous with-fixes design).
function band_curve(ref_id::AbstractString, S; mode::Symbol = :fixed)
    j = Dict("p10" => 1, "p50" => 2, "p90" => 3)[last(split(ref_id, '_'))]
    if mode === :monolix
        return S.times, S.median[:, j]                 # bin nodes joined by straight segments in (day, log10 y)
    elseif mode === :scheduled
        idx = [findfirst(==(t_), S.times) for t_ in SCHEDULED]
        return S.times[idx], S.median[idx, j]
    elseif startswith(ref_id, "fig3a_right")
        sel = findall(t_ -> t_ <= 30.0, S.times)
        return S.times[sel], S.median[sel, j]
    else
        idx = [findfirst(==(n), S.times) for n in NODES_LEFT]
        return S.times[idx], S.median[idx, j]          # compare_curve interpolates linearly in (day, log10 y)
    end
end

# dense left-panel curve (not graded; the sensitivity to the chord construction)
function dense_left_curve(ref_id::AbstractString, S)
    j = Dict("p10" => 1, "p50" => 2, "p90" => 3)[last(split(ref_id, '_'))]
    sel = findall(t_ -> t_ >= 31.0, S.times)
    return S.times[sel], S.median[sel, j]
end

function benchmark_setups()
    model = population_model()
    subj = cohort(3; times = TIMES_RIGHT)
    out = BenchmarkSetup[]
    # one typical individual (all etas zero) on the right-panel grid
    prob = QSPAtlas.base_problem(model, subj[1])
    push!(out, BenchmarkSetup(; name = "typical patient, Figure 3a right grid", prob = prob, saveat = TIMES_RIGHT))
    return out
end
