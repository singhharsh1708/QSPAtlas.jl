# Numerical convergence of every atlas score.
#
#   julia --startup-file=no --threads=4 --project=. benchmarks/convergence/run.jl [model_id ... | nostein | fornari_diag | report | check]
#
# `nostein` stands for every model except Stein (the CI job), followed by `report` and `check`. `check` alone rechecks the
# committed results. After running models the script checks them (see `check`) and exits non-zero if a curve flips
# pass/fail in any variant, a solve has a non-successful retcode (except the three documented Lindauer published bars),
# or, for a curve within 3 times its tolerance, an rms moves by more than 10 % of the tolerance.
#
# For each model the solves that feed the score are repeated with
#   own        the settings the model's protocol.jl uses (reproduces score.toml),
#   tol/100    abstol and reltol 100 times tighter,
#   altA, altB two other solvers (another family where the model allows it),
#   dense x4   4 times as many saved / observed points (Stein: 4 times as many observation times).
# Nothing in models/ or src/ is edited: the model modules are loaded as they are, their problem builders are
# reused, and the few lines that hard-code the solver are copied here with the solver as a parameter
# (each copy is marked "copy of" with the source function). Per curve and variant the script records rms,
# p95, max, pass/fail (the same compare_curve and tolerances as run_curves), the largest shift of the
# simulated curve against the "own" curve (in axis-box units, the unit of the score), and every solve's
# retcode and finiteness.
#
# Outputs (this folder): results/<id>_curves.csv, results/<id>_solves.csv, results/<id>_fidelity.csv,
# results.csv, solves.csv, fornari_human_reference.csv, tables.md (written by `report`).
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq
using Statistics, Printf, Random, TOML, DelimitedFiles
using SciMLBase: successful_retcode

const HERE = @__DIR__
const RES = joinpath(HERE, "results")
mkpath(RES)

# ---- configuration -------------------------------------------------------------------------------

struct Cfg
    name::String
    alg::Any
    atol::Float64
    rtol::Float64
    dens::Int
end

algname(a) = string(nameof(typeof(a)))

# own settings (copied from each model's solve call) and two alternative solvers
const OWN = Dict(
    "ouerdani_2015"   => (alg = Tsit5(),    atol = 1e-10, rtol = 1e-10, altA = Rodas5P(),     altB = Vern9()),
    "hutchinson_2016" => (alg = Tsit5(),    atol = 1e-9,  rtol = 1e-9,  altA = Rodas5P(),     altB = Vern9()),
    "lindauer_2017"   => (alg = Rodas5P(),  atol = 1e-14, rtol = 1e-9,  altA = FBDF(),        altB = Rosenbrock23()),
    "fornari_2019"    => (alg = Rodas5P(),  atol = 1e-9,  rtol = 1e-8,  altA = Vern9(),       altB = FBDF()),
    "stroh_2019"      => (alg = Rodas5P(),  atol = 1e-10, rtol = 1e-8,  altA = FBDF(),        altB = Rosenbrock23()),
    "stein_2019"      => (alg = Vern9(),    atol = 1e-10, rtol = 1e-10, altA = Rodas5P(),     altB = Tsit5()),
)

function variants(id)
    o = OWN[id]
    return [Cfg("own", o.alg, o.atol, o.rtol, 1),
        Cfg("tol/100", o.alg, o.atol / 100, o.rtol / 100, 1),
        Cfg("altA " * algname(o.altA), o.altA, o.atol, o.rtol, 1),
        Cfg("altB " * algname(o.altB), o.altB, o.atol, o.rtol, 1),
        Cfg("dense x4", o.alg, o.atol, o.rtol, 4)]
end

# ---- records -------------------------------------------------------------------------------------

const CURVES = NamedTuple[]
const SOLVES = NamedTuple[]

load_model(id) = (m = Module(Symbol("Conv_", id)); Base.include(m, joinpath(models_dir(), id, "reproduce.jl")); m)

fixkey(fixes) = join(sort([f.id for f in fixes]), ",")

# One solve, with retcode, timing and exception capture. `f()` returns the solution.
function timed_solve(id, run, key, cfg, f)
    t = @elapsed sol = try
        f()
    catch e
        e
    end
    if sol isa Exception
        push!(SOLVES, (; model = id, run, key, variant = cfg.name, retcode = "Exception: " * first(replace(sprint(showerror, sol), r"\s+" => " "), 120),
            ok = false, finite = false, npoints = 0, seconds = t))
        return nothing
    end
    ok = successful_retcode(sol)
    fin = try all(u -> all(isfinite, u), sol.u) catch; false end
    push!(SOLVES, (; model = id, run, key, variant = cfg.name, retcode = string(sol.retcode), ok, finite = fin, npoints = length(sol.t), seconds = t))
    return sol
end

# linear interpolation of (xv, yv) at x (duplicate x, as saved around a dose, keep the later value)
function lininterp(xv, yv, x)
    keep = [i for i in eachindex(xv) if i == length(xv) || xv[i + 1] != xv[i]]
    xs = xv[keep]; ys = yv[keep]
    return map(x) do xi
        j = searchsortedlast(xs, xi)
        j < 1 && return (xi >= xs[1] - 1e-9 ? ys[1] : NaN)
        j >= length(xs) && return (xi <= xs[end] + 1e-9 ? ys[end] : NaN)
        w = (xi - xs[j]) / (xs[j + 1] - xs[j])
        return ys[j] * (1 - w) + ys[j + 1] * w
    end
end

# largest change of the simulated curve against the own curve: in axis-box units (what the score sees) and relative
function curve_shift(ref, xb, yb, xv, yv)
    yvi = (length(xv) == length(xb) && all(xv .== xb)) ? yv : lininterp(xv, yv, xb)
    scale = maximum(abs, filter(isfinite, yb); init = 0.0)
    dbox = 0.0; drel = 0.0
    for i in eachindex(xb)
        a, b = yb[i], yvi[i]
        ua = isfinite(a) && QSPAtlas.usable(ref.yaxis, a); ub = isfinite(b) && QSPAtlas.usable(ref.yaxis, b)
        (ua || ub) || continue
        if ua != ub
            dbox = Inf; drel = Inf
            continue
        end
        dbox = max(dbox, abs(QSPAtlas.normalise(ref.yaxis, a) - QSPAtlas.normalise(ref.yaxis, b)))
        drel = max(drel, abs(a - b) / max(abs(a), 1e-9 * scale, floatmin()))
    end
    return dbox, drel
end

const OWNCURVE = Dict{Any, Any}()   # (id, run, curve) => (x, y, cmp)

function record_curve!(id, run, ref, cfg, tol, x, y, rc, finite)
    local cmp = nothing
    status = "ok"
    if x === nothing || !finite
        status = x === nothing ? "no result" : "non-finite output"
    end
    if x !== nothing
        cmp = try
            compare_curve(ref, x, y; tolerance = tol)
        catch e
            status = "compare failed: " * first(sprint(showerror, e), 80)
            nothing
        end
    end
    dbox = NaN; drel = NaN
    if cfg.name == "own"
        x !== nothing && (OWNCURVE[(id, run, ref.id)] = (x, y, cmp))
    elseif x !== nothing && haskey(OWNCURVE, (id, run, ref.id))
        xb, yb, _ = OWNCURVE[(id, run, ref.id)]
        dbox, drel = curve_shift(ref, xb, yb, x, y)
    end
    ownc = haskey(OWNCURVE, (id, run, ref.id)) ? OWNCURVE[(id, run, ref.id)][3] : nothing
    ownpass = ownc === nothing ? missing : ownc.passed
    margin = ownc === nothing ? NaN : max(ownc.rms / tol.rms, ownc.p95 / tol.p95)
    push!(CURVES, (; model = id, run, curve = ref.id, variant = cfg.name, status, retcode = rc, finite,
        rms = cmp === nothing ? NaN : cmp.rms, p95 = cmp === nothing ? NaN : cmp.p95, maxd = cmp === nothing ? NaN : cmp.max,
        passed = cmp === nothing ? false : cmp.passed, own_passed = ownpass, tol_rms = tol.rms, tol_p95 = tol.p95,
        own_margin = margin, shift_box = dbox, shift_rel = drel, n_points = x === nothing ? 0 : length(x)))
end

# ---- per-model simulation functions (solver as a parameter) --------------------------------------------
# Each returns (x, y, retcode::String, finite::Bool) for one reference curve, or throws NotRunnable.

# ouerdani_2015. copy of protocol.jl simulate_curve: Tsit5, abstol = reltol = 1e-10, saveat = 1000 points
function sim_ouerdani(M, id, run, ref, fixes, cfg)
    prob, mouse = M.problem_for(ref.id; fixes)
    n = (1000 - 1) * cfg.dens + 1
    ts = collect(range(0.0, prob.tspan[2]; length = n))
    sol = timed_solve(id, run, ref.id, cfg, () -> solve(prob, cfg.alg; abstol = cfg.atol, reltol = cfg.rtol, saveat = ts))
    sol === nothing && return (nothing, nothing, "Exception", false)
    y = try collect(Float64, endswith(ref.id, "_P") ? sol[M.TUMOUR.P] : sol[M.TUMOUR.K]) catch; fill(NaN, length(ts)) end
    length(y) == length(ts) || (y = vcat(y, fill(NaN, length(ts) - length(y))))
    return (mouse ? ts : ts ./ M.DAYS_PER_MONTH), y, string(sol.retcode), all(isfinite, y)
end

# hutchinson_2016. copy of protocol.jl simulate_scenario: Tsit5, 1e-9, breakpoint callback, dense output sampled on a 0.05 d grid
function sim_hutchinson(M, id, run, ref, fixes, cfg)
    p, out = M.scenario(ref.id)
    prob, cb = M.hutchinson_problem(p; fixes)
    kw = cb === nothing ? (;) : (; callback = cb)
    sol = timed_solve(id, run, ref.id, cfg, () -> solve(prob, cfg.alg; abstol = cfg.atol, reltol = cfg.rtol, kw...))
    sol === nothing && return (nothing, nothing, "Exception", false)
    grid = collect(M.T_START:(0.05 / cfg.dens):M.T_END)
    y = try collect(Float64, sol(grid; idxs = out === :T ? M.SYS.T : M.SYS.V).u) catch; fill(NaN, length(grid)) end
    return grid, y, string(sol.retcode), all(isfinite, y)
end

# lindauer_2017. copy of protocol.jl solve_human / log_ratio / cached_log_ratio / simulate_curve / simulate_bars
function lind_logratio(M, id, run, key, cfg, raw, mgkg, interval)
    prob, _ = M.human_problem(raw, mgkg, interval)
    sys = M.lindauer_system()
    kw = cfg.dens == 1 ? (; save_everystep = false, save_start = false) :
         (; saveat = collect(range(0.0, M.T_END_DAYS * 24; length = (366 - 1) * cfg.dens + 1)))
    sol = timed_solve(id, run, key, cfg, () -> solve(prob, cfg.alg; abstol = cfg.atol, reltol = cfg.rtol, kw..., dtmax = 24.0))
    sol === nothing && return NaN, "Exception"
    successful_retcode(sol) || return NaN, string(sol.retcode)
    v = sol(M.T_END_DAYS * 24; idxs = sys.TV)
    return (v > 0 ? log10(v / raw.W0) : NaN), string(sol.retcode)
end

const LCACHE = Dict{Any, Any}()
function lind_value(M, id, run, cfg, s, mgkg, name, mult)
    return get!(LCACHE, (id, run, cfg.name, s, mgkg, name, mult)) do
        raw = M.human_raw(; scenario = s.scenario, V_blood = s.V_blood, KdegPD1 = s.KdegPD1)
        if name != "base"
            f = M.TORNADO_FIELD[name]
            raw = merge(raw, NamedTuple{(f,)}((getfield(raw, f) * mult,)))
        end
        lind_logratio(M, id, run, string(mgkg, ":", name, ":", mult), cfg, raw, mgkg, s.interval_days)
    end
end

function sim_lindauer(M, id, run, ref, fixes, cfg)
    s = M.setup_from_fixes(fixes)
    if ref.id == "fig4_base_case"
        vs = [lind_value(M, id, run, cfg, s, M.DOSE_MGKG["5mgkg"], "base", 1.0), lind_value(M, id, run, cfg, s, M.DOSE_MGKG["0p5mgkg"], "base", 1.0)]
        return ref.x, first.(vs), join(unique(last.(vs)), "/"), all(isfinite, first.(vs))
    end
    m = match(r"^fig4_(5mgkg|0p5mgkg)_(x5|x0p2)$", ref.id)
    order = String.(split(ref.meta["order"], ", "))
    mult = m[2] == "x5" ? 5.0 : 0.2
    vs = [lind_value(M, id, run, cfg, s, M.DOSE_MGKG[m[1]], name, mult) for name in order]
    return ref.x, first.(vs), join(sort(unique(last.(vs))), "/"), all(isfinite, first.(vs))
end

# fornari_2019. copy of protocol.jl solved / case_solution / series: Rodas5P, abstol 1e-9, reltol 1e-8, saveat every 0.05 d
const FPROB = Dict{Any, Any}()
function sim_fornari(M, id, run, ref, fixes, cfg)
    case, v = M.TARGET_MAP[ref.id]
    dense = 0.05 / cfg.dens
    pkey = (case, fixkey(fixes), cfg.dens)
    s = get!(FPROB, pkey) do
        case === :rat ? M.rat_problem(fixes; dense) :
        case === :human_a ? M.human_problem(fixes; slopes = :rat, dense) : M.human_problem(fixes; slopes = :table2, dense)
    end
    skey = (:sol, id, run, case, cfg.name)
    sol = get!(FPROB, skey) do
        timed_solve(id, run, string(case), cfg, () -> solve(s.prob, cfg.alg; abstol = cfg.atol, reltol = cfg.rtol, saveat = s.saveat))
    end
    sol === nothing && return (nothing, nothing, "Exception", false)
    y = try [u for u in sol[getproperty(s.sys, v)]] catch; Float64[] end
    return sol.t ./ 24, y, string(sol.retcode), !isempty(y) && all(isfinite, y)
end

# stroh_2019. copy of protocol.jl run_spec / simulate_curve: Rodas5P, abstol 1e-10, reltol 1e-8, 6000 saved points per run
const SPROB = Dict{Any, Any}()
function sim_stroh(M, id, run, ref, fixes, cfg)
    sp, out = M.spec_for(ref.id)
    has(i) = any(f -> f.id == i, fixes)
    if startswith(String(out), "Upt_")
        has("readout_uptake_internalised") || throw(QSPAtlas.NotRunnable("Figure 4d,e 'Total Drug (nmol)' is not defined"))
    elseif has("unbound_readout") && out in (:Ctot_peripheral, :Ctot_tumor)
        out = Symbol(replace(String(out), "Ctot_" => "Cfree_"))
    else
        has("readout_total_all_forms") || throw(QSPAtlas.NotRunnable("'Total Drug' is not defined"))
    end
    r = get!(SPROB, (sp, fixkey(fixes), cfg.name, run)) do
        prob, cb, ts = M.make_problem(sp; fixes, npts = (6000 - 1) * cfg.dens + 1)
        sol = timed_solve(id, run, string(sp.species, ":", sp.molecule, ":", sp.dose), cfg,
            () -> solve(prob, cfg.alg; callback = cb, saveat = ts, abstol = cfg.atol, reltol = cfg.rtol))
        (; sol, ts)
    end
    r.sol === nothing && return (nothing, nothing, "Exception", false)
    y = try collect(Float64, r.sol(r.ts; idxs = getproperty(M.get_system(), out)).u) catch; fill(NaN, length(r.ts)) end
    return r.ts, y, string(r.sol.retcode), all(isfinite, y)
end

const SIMFN = Dict("ouerdani_2015" => sim_ouerdani, "hutchinson_2016" => sim_hutchinson, "lindauer_2017" => sim_lindauer,
    "fornari_2019" => sim_fornari, "stroh_2019" => sim_stroh)

# ---- generic driver ---------------------------------------------------------------------------------

function load_refs(id)
    return sort(collect(values(load_references(joinpath(models_dir(), id, "reference")))); by = r -> r.id)
end

# fidelity of the own variant against the committed score.toml and against the model's own simulate_curve
function fidelity(id, M, runs, refs)
    sc = read_score(joinpath(models_dir(), id, "score.toml"))
    rows = NamedTuple[]
    for (run, label, fixes) in runs
        tomlcurves = Dict(c["id"] => c for c in get(get(sc, run, Dict()), "curves", []))
        for ref in refs
            haskey(tomlcurves, ref.id) || continue
            key = (id, label, ref.id)
            haskey(OWNCURVE, key) || continue
            _, _, cmp = OWNCURVE[key]
            cmp === nothing && continue
            tc = tomlcurves[ref.id]
            push!(rows, (; model = id, run = label, curve = ref.id, rms = cmp.rms, toml_rms = tc["rms"], p95 = cmp.p95, toml_p95 = tc["p95"],
                passed = cmp.passed, toml_passed = tc["passed"],
                rms_rel_diff = abs(cmp.rms - tc["rms"]) / max(abs(tc["rms"]), 1e-12)))
        end
    end
    return rows
end

# a pair (score.toml key, label used here, fixes)
function model_runs(M)
    r = [("published", "published", M.PUBLISHED_FILLS)]
    isempty(M.FIXES) || push!(r, ("fixed", "fixed", M.FIXES))
    return r
end

function write_csv(path, rows)
    isempty(rows) && return
    ks = collect(keys(rows[1]))
    open(path, "w") do io
        println(io, join(ks, ","))
        for r in rows
            println(io, join([let v = r[k]; v isa AbstractString ? "\"" * replace(v, "\"" => "'") * "\"" : (v isa AbstractFloat ? @sprintf("%.10g", v) : string(v)) end for k in ks], ","))
        end
    end
end

function run_model(id)
    @info "model" id
    empty!(CURVES); empty!(SOLVES)
    M = load_model(id)
    id == "stein_2019" && return Base.invokelatest(run_stein, id, M)
    refs = load_refs(id)
    runs = model_runs(M)
    sim = SIMFN[id]
    for (key, label, fixes) in runs
        for cfg in variants(id)
            for ref in refs
                tol = get(M.TOLERANCES, ref.id, M.TOLERANCE)
                res = try
                    Base.invokelatest(sim, M, id, label, ref, fixes, cfg)
                catch e
                    e isa QSPAtlas.NotRunnable || rethrow()
                    nothing
                end
                res === nothing && continue          # NotRunnable: this run has no curves (as in run_curves)
                x, y, rc, fin = res
                record_curve!(id, label, ref, cfg, tol, x, y, rc, fin)
            end
            @info "done" id run = label variant = cfg.name solves = length(SOLVES)
        end
    end
    write_csv(joinpath(RES, "$(id)_curves.csv"), CURVES)
    write_csv(joinpath(RES, "$(id)_solves.csv"), SOLVES)
    write_csv(joinpath(RES, "$(id)_fidelity.csv"), fidelity(id, M, runs, refs))
    return nothing
end

# ---- stein_2019 (Monte Carlo) -------------------------------------------------------------------------
# The cohort simulation of protocol.jl, with three changes that make the variants comparable:
#  - solver and tolerance are parameters (population_model already takes both);
#  - the residual error is applied here from the SAME draws the production code uses (stream_rng(seed + block,
#    replicate, id, :residual), one randn per baseline observation time, in order), so that a variant with more
#    observation times sees identical residuals at the baseline times. The own variant is checked to reproduce
#    the production summary exactly;
#  - the percentile summary is the same type 7 quantile over the 90 patients, median over replicates.
const TIMES4 = sort(unique(vcat(collect(0.0625:0.0625:30.0), collect(31.0:0.25:190.0),
    [1.0, 2.0, 7.0, 11.0, 14.0, 21.0, 28.0, 3 * 365.25 / 12, 6 * 365.25 / 12])))   # the baseline TIMES, refined 4 times, plus the visits

function stein_summary(M, cfg; replicates, seed, nsub = M.N_PATIENTS)
    model = M.population_model(; solver = cfg.alg, tol = cfg.atol)
    tsim = cfg.dens == 1 ? M.TIMES : TIMES4
    subjects = M.cohort(nsub; times = tsim)
    bt = M.TIMES
    idx = [searchsortedfirst(tsim, t) for t in bt]
    all(i -> tsim[idx[i]] == bt[i], eachindex(bt)) || error("baseline times are not contained in the simulation times")
    nt = length(tsim); nb = length(bt)
    ids = [s.id for s in subjects]
    σ = M.RESIDUAL_SD_LOG
    Q = fill(NaN, nb, 3, replicates)
    nfail = 0; nnonfinite = 0; nrows = 0
    secs = @elapsed begin
        k0 = 0; b = 0
        while k0 < replicates
            m = min(M.CHUNK, replicates - k0)
            res = simulate(model, subjects; replicates = m, seed = seed + b, residual = false)
            ip = reshape([o.ipred for o in res.obs], nt, nsub, m)       # rows are ordered replicate, subject, time
            res.obs[1].time == tsim[1] && res.obs[nt].time == tsim[nt] && res.obs[nt + 1].time == tsim[1] || error("unexpected row order")
            nfail += count(s -> s.retcode !== :Success, res.subjects)
            nnonfinite += count(!isfinite, ip)
            nrows += length(ip)
            v = Vector{Float64}(undef, nsub)
            for r in 1:m
                eps = Matrix{Float64}(undef, nb, nsub)
                for j in 1:nsub
                    rng = QSPAtlas.stream_rng(seed + b, r, ids[j], :residual)
                    for i in 1:nb
                        eps[i, j] = randn(rng)
                    end
                end
                for i in 1:nb
                    for j in 1:nsub
                        v[j] = ip[idx[i], j, r] * exp(σ * eps[i, j])
                    end
                    q = quantile(v, M.QS)
                    for (kq, qv) in enumerate(q)
                        Q[i, kq, k0 + r] = qv
                    end
                end
            end
            k0 += m; b += 1
        end
    end
    S = M.summarise(bt, Q)
    return (; S, nfail, nnonfinite, nrows, seconds = secs)
end

# deterministic part: ipred of sampled individuals against the closed form of Eq. 1 (exact)
function stein_deterministic(M, cfg; n = 300, seed = 7)
    model = M.population_model(; solver = cfg.alg, tol = cfg.atol)
    tsim = cfg.dens == 1 ? M.TIMES : TIMES4
    subjects = M.cohort(n; times = tsim)
    worst = 0.0; worst_t = NaN; nbad = 0
    for η in (:zero, :sample)
        res = simulate(model, subjects; η, residual = false, seed)
        nbad += count(s -> s.retcode !== :Success, res.subjects)
        nt = length(tsim)
        for (k, sub) in enumerate(res.subjects)
            e = sub.η; θ = M.TYPICAL
            Tm = θ.Tmax * exp(e.Tmax)
            p = (foldx = θ.foldx * exp(e.foldx), Tmax = Tm, Cmax = θ.Cmax * exp(e.Cmax), α = θ.α * exp(e.α),
                FB = θ.FB * exp(e.FB), β = θ.β * exp(e.β), F1 = 1.0, F2 = 1.0, t1 = Tm, t2 = Tm)
            for i in 1:nt
                o = res.obs[(k - 1) * nt + i]
                cf = M.closed_form(o.time, p)
                rel = abs(o.ipred - cf) / cf
                isfinite(rel) || (rel = Inf)
                rel > worst && (worst = rel; worst_t = o.time)
            end
        end
    end
    return (; worst_rel = worst, at_time = worst_t, nfail = nbad)
end

function run_stein(id, M)
    refs = load_refs(id)
    reps = parse(Int, get(ENV, "CONV_STEIN_REPS", string(M.REPLICATES)))   # default: the protocol's 800
    vs = variants(id)
    extra = [Cfg("reps x2", OWN[id].alg, OWN[id].atol, OWN[id].rtol, 1)]
    detrows = NamedTuple[]
    mcrows = NamedTuple[]
    # production summary, for the fidelity check of the re-implementation
    prod = M.simulate_design(; replicates = reps)
    modes = [("published", :scheduled, M.PUBLISHED_FILLS), ("fixed", :fixed, M.FIXES)]
    own_S = nothing
    all_cfgs = vcat([(c, reps, M.SEED) for c in vs], [(extra[1], 2 * reps, M.SEED)],
        [(Cfg("seed+1000", OWN[id].alg, OWN[id].atol, OWN[id].rtol, 1), reps, M.SEED + 1000)])
    for (cfg, nrep, seed) in all_cfgs
        @info "stein variant" cfg.name nrep
        r = stein_summary(M, cfg; replicates = nrep, seed)
        cfg.name == "own" && (own_S = r.S;
            push!(mcrows, (; what = "own vs production summary (max abs diff of median)", value = maximum(abs, r.S.median .- prod.median))))
        push!(SOLVES, (; model = id, run = "mc", key = "cohort $(nrep) replicates x $(M.N_PATIENTS) patients", variant = cfg.name,
            retcode = r.nfail == 0 ? "Success (all $(r.nrows ÷ length(cfg.dens == 1 ? M.TIMES : TIMES4)) trajectories)" : "$(r.nfail) trajectories not Success",
            ok = r.nfail == 0, finite = r.nnonfinite == 0, npoints = r.nrows, seconds = r.seconds))
        for (label, mode, fixes) in modes
            for ref in refs
                tol = get(M.TOLERANCES, ref.id, M.TOLERANCE)
                x, y = M.band_curve(ref.id, r.S; mode)
                record_curve!(id, label, ref, cfg, tol, x, y, r.nfail == 0 ? "Success" : "failed", r.nnonfinite == 0)
            end
        end
        # shift of the compared medians against own in units of the Monte Carlo standard error
        if own_S !== nothing && cfg.name != "own"
            j = [1, 2, 3]
            for (label, mode, fixes) in modes, ref in refs
                jq = Dict("p10" => 1, "p50" => 2, "p90" => 3)[last(split(ref.id, '_'))]
                idxs = mode === :scheduled ? [findfirst(==(t_), r.S.times) for t_ in M.SCHEDULED] :
                       startswith(ref.id, "fig3a_left") ? [findfirst(==(n), r.S.times) for n in M.NODES_LEFT] :
                       findall(t_ -> t_ <= 30.0, r.S.times)
                dl = [abs(log10(r.S.median[i, jq]) - log10(own_S.median[i, jq])) for i in idxs]
                z = filter(isfinite, [abs(r.S.median[i, jq] - own_S.median[i, jq]) / own_S.mc_se[i, jq] for i in idxs])
                axd = log10(ref.yaxis.hi / ref.yaxis.lo)
                push!(mcrows, (; what = "median shift vs own: $label $(ref.id) $(cfg.name)",
                    value = maximum(dl) / axd, extra = isempty(z) ? NaN : maximum(z)))
            end
        end
        if cfg.name in ("own", "tol/100", vs[3].name, vs[4].name, "dense x4")
            d = stein_deterministic(M, cfg)
            push!(detrows, (; model = id, variant = cfg.name, max_rel_error_vs_closed_form = d.worst_rel, at_time = d.at_time, failed_trajectories = d.nfail))
        end
    end
    write_csv(joinpath(RES, "$(id)_curves.csv"), CURVES)
    write_csv(joinpath(RES, "$(id)_solves.csv"), SOLVES)
    write_csv(joinpath(RES, "$(id)_deterministic.csv"), detrows)
    write_csv(joinpath(RES, "$(id)_mc_shifts.csv"), [merge((; what = r.what, value = r.value), (extra = get(r, :extra, NaN),)) for r in mcrows])
    # fidelity: own against score.toml
    runs = [("published", "published", M.PUBLISHED_FILLS), ("fixed", "fixed", M.FIXES)]
    write_csv(joinpath(RES, "$(id)_fidelity.csv"), fidelity(id, M, runs, refs))
    return nothing
end

# ---- Fornari human Figure 4b reference solve --------------------------------------------------------------
# benchmarks/solvers/results/summary.md: "reference solve failed". Repeat what benchmarks/solvers/run.jl try_reference does
# (Vern9 with maxiters 1e5, then Rodas5P with maxiters 1e7, both at abstol = reltol = 1e-12, callback passed again although
# the problem already carries it) on the problem benchmark_setups() returns, then vary solver, tolerance, maxiters, callback.
fornari_diag() = Base.invokelatest(fornari_diag_run, load_model("fornari_2019"))

function fornari_diag_run(M)
    rows = NamedTuple[]
    function attempt(case, alg, atol, rtol, maxiters, dup, prob, cb, sys)
        kw = (dup && cb !== nothing) ? (; callback = cb) : (;)
        t = @elapsed sol = try
            solve(prob, alg; abstol = atol, reltol = rtol, maxiters, kw...)
        catch e
            e
        end
        if sol isa Exception
            push!(rows, (; case, solver = algname(alg), abstol = atol, reltol = rtol, maxiters, callback_twice = dup, retcode = "Exception: " * first(replace(sprint(showerror, sol), r"\s+" => " "), 100),
                steps = 0, seconds = t, final_Neut = NaN, min_Neut = NaN, min_Plt = NaN, min_state = NaN))
        else
            fin = !isempty(sol.u) && all(u -> all(isfinite, u), sol.u)
            push!(rows, (; case, solver = algname(alg), abstol = atol, reltol = rtol, maxiters, callback_twice = dup, retcode = string(sol.retcode),
                steps = length(sol.t), seconds = t, final_Neut = fin ? sol[sys.Neut][end] : NaN,
                min_Neut = fin ? minimum(sol[sys.Neut]) : NaN, min_Plt = fin ? minimum(sol[sys.Plt]) : NaN,
                min_state = fin ? minimum(minimum.(sol.u)) : NaN))
        end
    end
    h = Base.invokelatest(M.human_problem, M.FIXES; slopes = :table2)
    b = Base.invokelatest(M.benchmark_setups)[2]       # the BenchmarkSetup the benchmark script uses
    @assert b.prob.tspan == h.prob.tspan
    # 1. the benchmark's own try_reference sequence
    attempt("benchmark try_reference", Vern9(), 1e-12, 1e-12, 10^5, true, b.prob, b.callback, h.sys)
    attempt("benchmark try_reference", Rodas5P(), 1e-12, 1e-12, 10^7, true, b.prob, b.callback, h.sys)
    # 2. sweep: solver x tolerance, callback once (as in reproduce) and twice (as in the benchmark)
    for alg in (Tsit5(), Vern9(), Rodas5P(), FBDF()), tol in (1e-6, 1e-8, 1e-10, 1e-12, 1e-13), dup in (false, true)
        attempt("sweep", alg, tol, tol, 10^7, dup, h.prob, h.callback, h.sys)
    end
    # 3. what the benchmark's maxiters cap does to explicit solvers at 1e-12
    for alg in (Tsit5(), Vern9())
        attempt("maxiters 1e5", alg, 1e-12, 1e-12, 10^5, true, h.prob, h.callback, h.sys)
    end
    # 4. all other variants of the human problem the protocol defines (rat slopes 4a, Eq 10 slopes, single dose, bolus), both reference solvers
    fx = M.FIXES
    nofix(i) = filter(f -> f.id != i, fx)
    for (nm, f) in ("fig4a rat slopes" => () -> M.human_problem(fx; slopes = :rat), "eq10 slopes" => () -> M.human_problem(fx; slopes = :eq10),
            "single dose" => () -> M.human_problem(nofix("human_schedule"); slopes = :table2), "bolus instead of infusion" => () -> M.human_problem(fx; slopes = :table2, infusion = 0.0))
        hh = Base.invokelatest(f)
        for alg in (Vern9(), Rodas5P())
            attempt(nm, alg, 1e-12, 1e-12, alg isa Vern9 ? 10^5 : 10^7, true, hh.prob, hh.callback, hh.sys)
        end
    end
    # 5. baselines anywhere in the printed Table 1 ranges (Neut0 2-7, Mono0 0.2-10, Ret0 40-115, RBC0 4100-5900), corners and midpoints
    for N in (2.0, 4.5, 7.0), Mo in (0.2, 5.1, 10.0), R in (40.0, 77.5, 115.0), B in (4100.0, 5000.0, 5900.0)
        phys = M.human_physiology(; Neut0 = N, Mono0 = Mo, Ret0 = R, RBC0 = B)
        hh = Base.invokelatest(() -> M.human_problem(M.FIXES; slopes = :table2, phys))
        for alg in (Vern9(), Rodas5P())
            attempt("baseline range N=$N Mo=$Mo R=$R B=$B", alg, 1e-12, 1e-12, alg isa Vern9 ? 10^5 : 10^7, true, hh.prob, hh.callback, hh.sys)
        end
    end
    write_csv(joinpath(HERE, "fornari_human_reference.csv"), rows)
    bad = filter(r -> !startswith(r.retcode, "Success"), rows)
    @info "fornari human reference" attempts = length(rows) not_success = length(bad)
    for r in bad
        @info "not successful" r.case r.solver r.abstol r.retcode
    end
    return rows
end

# ---- report ---------------------------------------------------------------------------------------------

function read_csv(path)
    isfile(path) || return NamedTuple[]
    data, hdr = readdlm(path, ','; header = true, quotes = true)
    return [NamedTuple{Tuple(Symbol.(hdr[:]))}(Tuple(data[i, :])) for i in 1:size(data, 1)]
end

mx(v) = isempty(v) ? NaN : maximum(v)
fmt(x) = x isa Real ? (isnan(x) ? "nan" : isinf(x) ? "inf" : @sprintf("%.3g", x)) : string(x)

function report()
    ids = collect(keys(OWN))
    sort!(ids)
    allc = NamedTuple[]; alls = NamedTuple[]
    for id in ids
        append!(allc, read_csv(joinpath(RES, "$(id)_curves.csv")))
        append!(alls, read_csv(joinpath(RES, "$(id)_solves.csv")))
    end
    write_csv(joinpath(HERE, "results.csv"), allc)
    write_csv(joinpath(HERE, "solves.csv"), alls)
    io = IOBuffer()
    println(io, "### Per model and variant: largest score change against the own settings\n")
    println(io, "`d rms` and `d p95` are the largest absolute change of rms / p95 over the curves of the run (axis-box units; the pass rule is rms <= 0.015 and p95 <= 0.03 unless a curve has its own tolerance). ",
        "`shift` is the largest change of any simulated point in axis-box units. `flips` counts curves whose pass/fail differs from own.\n")
    println(io, "| model | run | variant | curves | max d rms | max d p95 | max shift (box) | flips | failed solves | non-finite |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|")
    for id in ids, run in ("published", "fixed")
        cs = filter(r -> r.model == id && r.run == run, allc)
        isempty(cs) && continue
        own = Dict(r.curve => r for r in cs if r.variant == "own")
        for var in unique(r.variant for r in cs)
            var == "own" && continue
            vc = filter(r -> r.variant == var, cs)
            drms = [abs(r.rms - own[r.curve].rms) for r in vc if haskey(own, r.curve) && isfinite(r.rms) && isfinite(own[r.curve].rms)]
            dp95 = [abs(r.p95 - own[r.curve].p95) for r in vc if haskey(own, r.curve) && isfinite(r.p95) && isfinite(own[r.curve].p95)]
            sh = [r.shift_box for r in vc if r.shift_box isa Real && !isnan(r.shift_box)]
            flips = count(r -> haskey(own, r.curve) && r.passed != own[r.curve].passed, vc)
            sol = filter(r -> r.model == id && r.run in (run, "mc") && r.variant == var, alls)
            nfail = count(r -> !(r.ok === true || r.ok == "true"), sol)
            nnf = count(r -> !(r.finite === true || r.finite == "true"), sol)
            println(io, "| ", id, " | ", run, " | ", var, " | ", length(vc), " | ", fmt(mx(drms)), " | ", fmt(mx(dp95)), " | ",
                fmt(mx(sh)), " | ", flips, " | ", nfail, " | ", nnf, " |")
        end
    end
    println(io, "\n### Per curve: the three curves with the largest change of rms in each model and run\n")
    println(io, "Largest change over all variants; the full per-curve, per-variant table is results.csv.\n")
    println(io, "| model | run | curve | own rms | own p95 | max d rms | max d p95 | max shift (box) | variant of largest d rms |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for id in ids, run in ("published", "fixed")
        cs = filter(r -> r.model == id && r.run == run, allc)
        isempty(cs) && continue
        stats = NamedTuple[]
        for cv in unique(r.curve for r in cs)
            rs = filter(r -> r.curve == cv, cs)
            o = filter(r -> r.variant == "own", rs)
            isempty(o) && continue
            o = o[1]
            vs = [r for r in rs if r.variant != "own" && isfinite(r.rms) && isfinite(o.rms)]
            isempty(vs) && continue
            d = [abs(r.rms - o.rms) for r in vs]
            k = argmax(d)
            push!(stats, (; cv, o, drms = d[k], dp95 = mx([abs(r.p95 - o.p95) for r in vs]), shift = mx([r.shift_box for r in vs if r.shift_box isa Real && !isnan(r.shift_box)]), var = vs[k].variant))
        end
        sort!(stats; by = x -> -x.drms)
        for x in first(stats, 3)
            println(io, "| ", id, " | ", run, " | ", x.cv, " | ", fmt(x.o.rms), " | ", fmt(x.o.p95), " | ", fmt(x.drms), " | ", fmt(x.dp95), " | ", fmt(x.shift), " | ", x.var, " |")
        end
    end
    println(io, "\n### Curves whose pass/fail flips in some variant, or that sit within 1.5x of a threshold\n")
    println(io, "| model | run | curve | own rms | own p95 | own pass | own margin (max of rms/tol, p95/tol) | variants that flip |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for id in ids, run in ("published", "fixed")
        cs = filter(r -> r.model == id && r.run == run, allc)
        for cv in unique(r.curve for r in cs)
            rs = filter(r -> r.curve == cv, cs)
            o = filter(r -> r.variant == "own", rs)
            isempty(o) && continue
            o = o[1]
            fl = [r.variant for r in rs if r.variant != "own" && r.passed != o.passed]
            (isempty(fl) && !(0.67 <= o.own_margin <= 1.5)) && continue
            println(io, "| ", id, " | ", run, " | ", cv, " | ", fmt(o.rms), " | ", fmt(o.p95), " | ", o.passed, " | ", fmt(o.own_margin), " | ", isempty(fl) ? "none" : join(fl, "; "), " |")
        end
    end
    println(io, "\n### Solves that did not finish with a successful retcode or gave non-finite output\n")
    bad = filter(r -> !(r.ok === true || r.ok == "true") || !(r.finite === true || r.finite == "true"), alls)
    if isempty(bad)
        println(io, "None.")
    else
        println(io, "| model | run | solve | variant | retcode | finite |")
        println(io, "|---|---|---|---|---|---|")
        for r in bad
            println(io, "| ", r.model, " | ", r.run, " | ", r.key, " | ", r.variant, " | ", r.retcode, " | ", r.finite, " |")
        end
    end
    println(io, "\n### Reproduction of score.toml by the own variant\n")
    println(io, "| model | curves compared | max relative rms difference | pass/fail differences |")
    println(io, "|---|---|---|---|")
    for id in ids
        f = read_csv(joinpath(RES, "$(id)_fidelity.csv"))
        isempty(f) && continue
        println(io, "| ", id, " | ", length(f), " | ", fmt(maximum(r.rms_rel_diff for r in f)), " | ", count(r -> r.passed != r.toml_passed, f), " |")
    end
    write(joinpath(HERE, "tables.md"), String(take!(io)))
    return nothing
end

# ---- check (exit status) -------------------------------------------------------------------------------

# Solves that are known not to finish in the published Lindauer run: the simulated tumour falls to about 1e-11 mm^3 and
# every concentration scales with 1/TV (see models/lindauer_2017/report.md). They are scored as complete misses.
const ALLOWED_FAILED_SOLVES = Set([("lindauer_2017", "published", "0.5:gamma:5.0"), ("lindauer_2017", "published", "5.0:gamma:5.0"),
    ("lindauer_2017", "published", "5.0:SLtg:5.0")])
const RMS_MOVE_FRACTION = 0.1   # of the curve's rms tolerance
const MC_NOISE_VARIANTS = ("reps x2", "seed+1000")

istrue(x) = x === true || x == "true"

# Problems found in the results of the models `ids`: flips of pass/fail, failed solves, rms moves above 10 % of the
# tolerance for curves within NEAR_THRESHOLD times the tolerance (clearly failing curves cannot flip).
const NEAR_THRESHOLD = 3.0
function check_problems(ids)
    problems = String[]
    for id in ids
        curves = read_csv(joinpath(RES, "$(id)_curves.csv"))
        solves = read_csv(joinpath(RES, "$(id)_solves.csv"))
        isempty(curves) && push!(problems, "$id: no curves in results/$(id)_curves.csv")
        for r in solves
            istrue(r.ok) && continue
            (r.model, r.run, r.key) in ALLOWED_FAILED_SOLVES && continue
            push!(problems, "$id: solve $(r.run) $(r.key) [$(r.variant)] ended with $(r.retcode)")
        end
        for run in unique(r.run for r in curves), cv in unique(r.curve for r in curves if r.run == run)
            rs = filter(r -> r.run == run && r.curve == cv, curves)
            o = filter(r -> r.variant == "own", rs)
            isempty(o) && continue
            o = o[1]
            for r in rs
                r.variant == "own" && continue
                r.passed != o.passed && push!(problems, "$id: $run $cv flips pass/fail under $(r.variant) (own $(o.passed))")
                r.variant in MC_NOISE_VARIANTS && continue   # Stein's Monte Carlo noise references: they size the noise, not the solver error
                if isfinite(r.rms) && isfinite(o.rms)
                    # The rms-move limit protects results near the pass/fail line. A curve that fails
                    # by more than NEAR_THRESHOLD times its tolerance in both variants cannot flip, and
                    # its rms can legitimately depend on the output grid (for example where a
                    # simulation dips to zero on a log axis and the gap boundary moves with sampling).
                    near = min(r.rms, o.rms) <= NEAR_THRESHOLD * o.tol_rms
                    near && abs(r.rms - o.rms) > RMS_MOVE_FRACTION * o.tol_rms &&
                        push!(problems, "$id: $run $cv rms moves by $(fmt(abs(r.rms - o.rms))) under $(r.variant) (limit $(fmt(RMS_MOVE_FRACTION * o.tol_rms)))")
                elseif isfinite(r.rms) != isfinite(o.rms)
                    push!(problems, "$id: $run $cv rms is finite in only one of own and $(r.variant)")
                end
            end
        end
    end
    return problems
end

function check(ids)
    problems = check_problems(ids)
    if isempty(problems)
        println("convergence check passed for: ", join(ids, ", "))
        return nothing
    end
    println(stderr, "convergence check FAILED (", length(problems), " problems):")
    foreach(p -> println(stderr, "  ", p), problems)
    exit(1)
end

# ---- main ---------------------------------------------------------------------------------------------

function main(args)
    args = String[args...]
    if "nostein" in args
        filter!(!=("nostein"), args)
        args = vcat(sort(filter(!=("stein_2019"), collect(keys(OWN)))), args, ["report", "check"])
    end
    ids = isempty(args) ? sort(collect(keys(OWN))) : args
    ran = String[]
    for a in ids
        a in ("report", "fornari_diag", "check") && continue
        a in keys(OWN) || error("unknown model $a")
        run_model(a)
        push!(ran, a)
    end
    ("fornari_diag" in args || isempty(args)) && fornari_diag()
    ("report" in args || isempty(args)) && report()
    if !isempty(ran) || "check" in args
        check(isempty(ran) ? sort(collect(keys(OWN))) : ran)
    end
end

main(ARGS)
