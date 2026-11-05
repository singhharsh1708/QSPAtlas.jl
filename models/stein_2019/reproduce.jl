using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
import CairoMakie
using Random: Xoshiro

Base.include(@__MODULE__, joinpath(@__DIR__, "model.jl"))      # Base.include: a bare Module() has no include
Base.include(@__MODULE__, joinpath(@__DIR__, "protocol.jl"))

const ID = "stein_2019"
const CITATION = "Stein AM, et al. CPT Pharmacometrics Syst Pharmacol 2019;8:285-295"
const DOI = "10.1002/psp4.12388"
const PMCID = "PMC6539725"

# The paper's supplementary data set is not distributed with the repository. Optional local file, read from
# $QSPATLAS_PAPERS/stein_2019/supplement/ (default: a `papers` folder next to the repository, see tools/digitise/README.md).
const S1_CSV = joinpath(abspath(expanduser(get(ENV, "QSPATLAS_PAPERS", joinpath(@__DIR__, "..", "..", "..", "papers")))),
    "stein_2019", "supplement", "PSP4-8-285-s002.csv")

# Model check: the ODE (effector/memory compartments with ifelse switches) against the closed form of
# Eq. 1 / Eq. 2 / supplement Monolix code, for random individuals, with and without comedication.
function check_model(; n = 30, tol = 1.0e-10)
    worst = 0.0
    for comed in (false, true)
        m = population_model(; comedication = comed, tol)
        rng = Xoshiro(1)
        for trial in 1:n
            toci = comed && rand(rng) < 0.5 ? 1 + 20 * rand(rng) : Inf
            ster = comed && rand(rng) < 0.5 ? 0.1 + 20 * rand(rng) : Inf
            s = Subject("a"; regimen = Regimen(Dose[]), obstimes = collect(0.25:0.25:60.0), covariates = (toci = toci, ster = ster))
            r = simulate(m, s; rng = Xoshiro(trial), residual = false)
            η = r.subjects[1].η
            θ = TYPICAL
            Tm = θ.Tmax * exp(η.Tmax)
            t1, f1, t2, f2 = !comed ? (Tm, 1.0, Tm, 1.0) :
                toci <= ster ? (min(toci, Tm), θ.Ftoci, min(ster, Tm), θ.Fster) : (min(ster, Tm), θ.Fster, min(toci, Tm), θ.Ftoci)
            p = (foldx = θ.foldx * exp(η.foldx), Tmax = Tm, Cmax = θ.Cmax * exp(η.Cmax), α = θ.α * exp(η.α),
                FB = θ.FB * exp(η.FB), β = θ.β * exp(η.β), F1 = f1, F2 = f2, t1 = t1, t2 = t2)
            for o in r.obs
                worst = max(worst, abs(o.ipred - closed_form(o.time, p)) / closed_form(o.time, p))
            end
        end
    end
    return worst
end

# Simulation summary of the graded design (base model, residual error included, no BLQ censoring).
function simulate_design(; replicates = REPLICATES, blq = nothing, model = population_model(), subjects = cohort(), seed = SEED)
    times, Q = replicate_percentiles(model, subjects; replicates, chunk = CHUNK, seed, blq)
    return summarise(times, Q)
end

# Monte Carlo standard error of the compared median, in decades and in axis-box units, over the
# scored reference x of each band.
function mc_report(S, refs; mode = :fixed)
    out = Dict{String, NamedTuple}()
    for ref in refs
        j = Dict("p10" => 1, "p50" => 2, "p90" => 3)[last(split(ref.id, '_'))]
        left = startswith(ref.id, "fig3a_left")
        idx = mode === :scheduled ? [findfirst(==(t_), S.times) for t_ in SCHEDULED] :
            left ? [findfirst(==(n), S.times) for n in NODES_LEFT] :
            [i for (i, t_) in enumerate(S.times) if t_ <= 30.0 && any(abs(t_ - x) <= 0.13 for x in ref.x)]
        dec = [S.mc_se[i, j] / (S.median[i, j] * log(10)) for i in idx]
        axis_decades = log10(ref.yaxis.hi / ref.yaxis.lo)
        out[ref.id] = (max_decades = maximum(dec), max_box = maximum(dec) / axis_decades, mean_box = sum(dec) / length(dec) / axis_decades)
    end
    return out
end

# Consistency of the chunked percentile code with QSPAtlas `vpc` (Part 4) on a small case: same cohort,
# same seed, one chunk. Returns the largest absolute difference of the median (should be 0).
function check_vpc_consistency(; replicates = 20, times = [2.0, 8.0, 12.0, 28.0, 100.0])
    model = population_model()
    subj = cohort(N_PATIENTS; times)
    v = vpc(model, subj; endpoint = :C, bins = :obstimes, replicates, seed = SEED, pi_level = PI_LEVEL)
    S = summarise(replicate_percentiles(model, subj; replicates, chunk = replicates, seed = SEED)...)
    worst = 0.0
    for r in v.rows
        i = findfirst(==(r.t_mid), S.times); j = findfirst(==(r.q), collect(QS))
        worst = max(worst, abs(r.median - S.median[i, j]), abs(r.lo - S.lo[i, j]), abs(r.hi - S.hi[i, j]), abs(r.mc_se - S.mc_se[i, j]))
    end
    return worst
end

# Compare a summary with the references (same comparison as run_curves, without a Run).
function compare_summary(S, refs; left = :chord, mode = :fixed, tolerance = TOLERANCE)
    return [compare_curve(ref, (left === :dense && startswith(ref.id, "fig3a_left") ? dense_left_curve(ref.id, S) : band_curve(ref.id, S; mode))...; tolerance) for ref in refs]
end

# ---- sensitivity runs (not graded; reported in report.md) ---------------------------------------
# Comedication design is ASSUMED: each patient gets tocilizumab with probability 0.36 and corticosteroids
# with probability 0.26 (Table 2), independently, first-dose day log-normal with the Table 2 median
# (5.7 d / 7.5 d) truncated to the Table 2 range (1-27 d / 0.11-170 d). The paper gives no joint design.
function comedication_sampler(rng)
    draw(median_, σ, lo, hi) = clamp(median_ * exp(σ * randn(rng)), lo, hi)
    toci = rand(rng) < 0.36 ? draw(5.7, 0.7, 1.0, 27.0) : Inf
    ster = rand(rng) < 0.26 ? draw(7.5, 1.2, 0.11, 170.0) : Inf
    return (toci = toci, ster = ster)
end

function sensitivity_runs(refs; replicates = REPLICATES)
    runs = Pair{String, Any}[]
    push!(runs, "omega read as variance (sd = sqrt)" => simulate_design(; replicates, model = population_model(omega = :variance)))
    push!(runs, "BLQ censored to LLOQ 50" => simulate_design(; replicates, blq = (lloq = LLOQ, handling = :censor_to_lloq)))
    push!(runs, "comedication model (assumed timing), Ftoci 1.2, Fster 1.0" =>
        simulate_design(; replicates, model = population_model(comedication = true), subjects = cohort(; comedication = comedication_sampler)))
    push!(runs, "residual sd halved (0.28 instead of 0.56)" => simulate_design(; replicates, model = population_model(sigma = 0.28)))
    return runs
end

# ---- figures ------------------------------------------------------------------------------
# Figure 3a in its layout: simulated percentile intervals (blue 10th/90th, pink 50th) with the median
# of each simulated percentile; one simulated cohort at the nominal visit days as dots, simulated BLQ as
# red asterisks. The observed data and the empirical percentile lines of the paper are not shown.
function plot_fig3(path, S)
    M = CairoMakie
    model = population_model()
    nominal = [1.0, 2.0, 7.0, 11.0, 14.0, 21.0, 28.0, 91.0, 182.0, 273.0, 365.0]
    one = simulate(model, cohort(N_PATIENTS; times = nominal); seed = 7)
    fig = M.Figure(size = (1100, 420))
    for (j, (xlim, ylim, xl)) in enumerate((((0, 380), (20, 4.0e5), "Days"), ((0, 30), (20, 4.0e5), "Days")))
        ax = M.Axis(fig[1, j]; xlabel = xl, ylabel = j == 1 ? "Tisagenlecleucel transgene (copies/ug)" : "", yscale = log10, limits = (xlim, ylim),
            title = j == 1 ? "Simulated VPC bands (this reproduction)" : "zoom, days 0-30")
        sel = findall(t_ -> t_ <= xlim[2], S.times)
        for (k, col) in ((1, :steelblue), (3, :steelblue), (2, :hotpink))
            M.band!(ax, S.times[sel], S.lo[sel, k], S.hi[sel, k]; color = (col, 0.35))
            M.lines!(ax, S.times[sel], S.median[sel, k]; color = col, linewidth = 2)
        end
        pts = [(o.time, o.dv) for o in one.obs if o.time <= xlim[2]]
        ok = [p for p in pts if p[2] >= LLOQ]
        bq = [p for p in pts if p[2] < LLOQ]
        M.scatter!(ax, first.(ok), last.(ok); color = (:navy, 0.5), markersize = 4)
        isempty(bq) || M.scatter!(ax, first.(bq), last.(bq); color = :crimson, marker = :star5, markersize = 6)
        M.hlines!(ax, [LLOQ]; color = :gray, linestyle = :dash)
    end
    mkpath(dirname(path))
    M.save(path, fig)
    return path
end

# Local diagnostic (skipped when the file is absent, e.g. in CI): empirical percentiles of the paper's
# Data set S1 (a simulated data set with the real design) against the simulated bands. Nothing is
# graded or committed from it; only per-bin percentiles are used, never the rows.
function s1_diagnostic(S, path; csv = S1_CSV)
    isfile(csv) || return nothing
    lines = readlines(csv)[2:end]
    rows = [split(replace(l, '"' => ""), ',') for l in lines if !isempty(strip(l))]
    tt = [parse(Float64, r[2]) for r in rows]
    yy = [exp(parse(Float64, r[3])) for r in rows]
    edges = [0.0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 14, 16, 18, 21, 24, 28, 35, 60, 100, 190, 400]
    M = CairoMakie
    fig = M.Figure(size = (700, 460))
    ax = M.Axis(fig[1, 1]; xlabel = "Days", ylabel = "copies/ug", yscale = log10, limits = (0, 60, 5, 1.0e7),
        title = "local only: Data set S1 percentiles (lines) vs simulated bands")
    for (k, col) in ((1, :steelblue), (3, :steelblue), (2, :hotpink))
        sel = findall(t_ -> t_ <= 60, S.times)
        M.band!(ax, S.times[sel], S.lo[sel, k], S.hi[sel, k]; color = (col, 0.35))
        qv = Float64[]; xm = Float64[]
        for b in 1:(length(edges) - 1)
            v = [y for (t_, y) in zip(tt, yy) if edges[b] <= t_ < edges[b + 1]]
            length(v) >= 5 || continue
            push!(xm, sum(t_ for t_ in tt if edges[b] <= t_ < edges[b + 1]) / count(t_ -> edges[b] <= t_ < edges[b + 1], tt))
            push!(qv, quantile(v, (0.1, 0.5, 0.9)[k]))
        end
        M.lines!(ax, xm, qv; color = :black, linewidth = 1.5)
    end
    mkpath(dirname(path))
    M.save(path, fig)
    return path
end

function reproduce(; plots = true, figdir = joinpath(@__DIR__, "figures"), replicates = REPLICATES, sensitivity = false)
    refs = sort(collect(values(load_references(joinpath(@__DIR__, "reference")))); by = r -> r.id)
    S = simulate_design(; replicates)
    published = run_curves(ref -> band_curve(ref.id, S; mode = band_mode(PUBLISHED_FILLS)), refs, "as published";
        fixes = PUBLISHED_FILLS, tolerance = TOLERANCE, tolerances = TOLERANCES)
    fixed = run_curves(ref -> band_curve(ref.id, S; mode = band_mode(FIXES)), refs, "with fixes"; fixes = FIXES, tolerance = TOLERANCE, tolerances = TOLERANCES)
    rep = Reproduction(; id = ID, citation = CITATION, doi = DOI, pmcid = PMCID, published, fixed)
    write_score(joinpath(@__DIR__, "score.toml"), rep)
    global LAST = (; S, refs, mc = mc_report(S, refs), mc_sched = mc_report(S, refs; mode = :scheduled))
    if sensitivity
        global SENS = [name => (compare_summary(S2, refs; mode = :scheduled), compare_summary(S2, refs), compare_summary(S2, refs; left = :dense)) for (name, S2) in sensitivity_runs(refs; replicates)]
    end
    if plots
        mkpath(figdir)
        plot_reproduction(rep, joinpath(figdir, "comparison.png"); ncols = 3, run = :published)
        plot_reproduction(rep, joinpath(figdir, "comparison_fixed.png"); ncols = 3, run = :fixed)
        plot_fig3(joinpath(figdir, "fig3.png"), S)
        s1_diagnostic(S, joinpath(figdir, "digitise_check_s1_diagnostic.png"))
    end
    return rep
end

(abspath(PROGRAM_FILE) == (@__FILE__)) && display(reproduce())
