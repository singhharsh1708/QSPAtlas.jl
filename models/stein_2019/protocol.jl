# Protocol, targets, tolerances, fixes for stein_2019 (Figure 3a, simulated percentile bands).
# TARGETS AND DESIGN FROZEN 2026-11-04 14:50 IST, before any simulation was compared with the
# bands (reference/*.csv, digitise.toml, the design constants and the tolerance below).
# No simulation had been run against the targets at that time (only the closed-form check of the
# model equations, which does not touch the bands). REVISED after review (see the design section and
# report.md): the as-published run now uses the scheduled visits; targets and tolerance unchanged except
# two misassigned points of fig3a_right_p90 removed by a documented rule.
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
#   "with fixes"   = the same simulation read on a denser design: a dense grid (0.25 d on (0, 30])
#       for the right panel and the left-panel bin nodes of the band (28, 68, 190 d) joined by
#       chords. See FIXES for each choice and its evidence.
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

# ---- fixes ----------------------------------------------------------------------------------
# Every design choice that is not stated by the paper. `evidence` says what supports it. The first
# five are used by both runs (the as-published run is not definable without them); the last two
# are what distinguishes the "with fixes" run from "as published".
const FIXES = [
    Fix("all_patients_all_times", :missing,
        "Every virtual patient is observed at every time of the design (no missing visits, no unscheduled samples).";
        evidence = "Methods 'Sample analysis' lists scheduled visits (days 1, 2, 7, 11, 14, 21, 28, then 3-monthly) and Results 836 measurements in 90 patients (9.3 per patient, 43 BLQ); who was sampled when is not given. Used by both runs."),
    Fix("residual_included", :ambiguity,
        "Percentiles are taken over observations WITH residual error (log-normal, sd 0.56).";
        evidence = "Fig 3a plots observations against simulated percentiles of observations (a VPC); Table 1 residual error a = 0.56. The caption does not say whether ipred or dv percentiles are shown. Used by both runs."),
    Fix("blq_kept", :ambiguity,
        "Simulated values below the LLOQ (50 copies/ug) are kept, not replaced by the LLOQ.";
        evidence = "Caption: red asterisks are 'simulated data points for the measurements that were BLQ'. This shows simulated values are plotted; it does not say how they enter the percentiles. Censoring to the LLOQ is a sensitivity run. Used by both runs."),
    Fix("quantile_type7", :missing,
        "Percentiles are Hyndman-Fan type 7 (linear interpolation of order statistics) over the 90 patients, median over replicates.";
        evidence = "Quantile convention of the VPC is not stated; type 7 is the R/Monolix-style default and is the QSPAtlas `percentiles` default. Used by both runs."),
    Fix("base_model_no_comedication", :ambiguity,
        "Base model (Eq. 1) without covariates and without the tocilizumab/corticosteroid effect (F = 1).";
        evidence = "Results 'Analytical model': 'The base model characterized the data well, as shown by the visual predictive check and individual fits in Figure 3'; Results 'Covariate analysis': the final model is the base model. Methods 'Covariate model' says the full covariate model is the final model, and Table 1 lists Ftoci = 1.2, so the VPC may have used comedication regressors; their per-patient values are not published. Comedication is a sensitivity run with assumed timings. Used by both runs."),
    Fix("dense_grid", :missing,
        "Right panel read on a dense grid (0.25 d) instead of only at the scheduled visits.";
        evidence = "Monolix VPC bins come from the data's actual (unscheduled, irregular) sample times, which are not published; the band between visits is drawn from those bins. A dense grid removes the coarse chord between scheduled visits, which is not how the band was built. Only in the 'with fixes' run."),
    Fix("left_bin_nodes", :from_target,
        "Left panel simulated percentile taken at three bin nodes (28, 68, 190 d) and joined by chords.";
        evidence = "Nodes are read from the kinks of the authors' band edges in Figure 3a (digitise.toml); the bins themselves are not stated, so this input comes from the target figure (circular, the dagger case). Only in the 'with fixes' run."),
]

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
        theta = TYPICAL, omegas = OMEGA, sigma = RESIDUAL_SD_LOG)
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
        observed = (C = CART.C,), error = (C = ExponentialError(sigma),), tspan = TSPAN,
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

# curve of one band target from a summary. mode :scheduled = scheduled visits joined by chords (as
# published); :fixed = right panel dense, left panel bin nodes joined by chords (with fixes).
function band_curve(ref_id::AbstractString, S; mode::Symbol = :fixed)
    j = Dict("p10" => 1, "p50" => 2, "p90" => 3)[last(split(ref_id, '_'))]
    if mode === :scheduled
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
