# Solver work-precision benchmarks on the exact problems used to reproduce each paper.
#
#   julia --project=benchmarks/solvers benchmarks/solvers/run.jl [model_id ...]
#
# For every BenchmarkSetup a model's protocol.jl exposes, a reference solution is
# computed with Vern9 (non-stiff) or Rodas5P at 1e-12 (falls back to the other when
# the first fails), then each solver is timed over a tolerance sweep. Error is the
# dense L2 error over the time span (interpolants compared; see benchmark_problem for why
# not saved points). Results: results/<model>__<setup>.png,
# results/summary.md (fastest solver reaching dense L2 error <= 1e-6 x RMS of the solution) and
# results/summary.toml.
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq, DiffEqDevTools
using OrdinaryDiffEqSDIRK: KenCarp4, TRBDF2
using OrdinaryDiffEqLowOrderRK: BS3
using OrdinaryDiffEqBDF: FBDF, QNDF
using OrdinaryDiffEqRosenbrock: Rodas5P, Rodas4P, Rosenbrock23
using OrdinaryDiffEqVerner: Vern7, Vern9
using OrdinaryDiffEqTsit5: Tsit5
using CairoMakie, TOML, Printf, Statistics
using SciMLBase: successful_retcode

const OUT = joinpath(@__DIR__, "results")
mkpath(OUT)

const SOLVERS = [
    Dict(:alg => Tsit5()), Dict(:alg => Vern7()), Dict(:alg => BS3()),
    Dict(:alg => Rodas5P()), Dict(:alg => Rodas4P()), Dict(:alg => Rosenbrock23()),
    Dict(:alg => FBDF()), Dict(:alg => QNDF()), Dict(:alg => KenCarp4()), Dict(:alg => TRBDF2()),
]
const ABSTOLS = 1.0 ./ 10.0 .^ (6:11)
const LEVEL_SECONDS = 10.0  # wall time allowed for one tolerance level (3 timed solves plus the error estimate)
const PROBE_SECONDS = 5.0   # wall time of one solve at the loosest tolerance, after compilation
const RELTOLS = 1.0 ./ 10.0 .^ (3:8)

function try_reference(b)
    for alg in (Vern9(), Rodas5P())
        kw = b.callback === nothing ? (;) : (; callback = b.callback)
        sol = try
            solve(b.prob, alg; abstol = 1e-12, reltol = 1e-12, maxiters = 10^7, kw...)
        catch
            nothing
        end
        sol !== nothing && successful_retcode(sol) && return alg, sol
    end
    return nothing, nothing
end

function bench_one(id, b)
    refalg, refsol = try_reference(b)
    refalg === nothing && return (; id, name = b.name, status = "reference solve failed", rows = [], scale = NaN)
    # typical magnitude of the saved solution, so the error target is relative
    scale = sqrt(mean(abs2, reduce(vcat, refsol(range(refsol.t[1], refsol.t[end]; length = 2000)).u)))
    rows = []
    fig = Figure(size = (720, 480))
    ax = CairoMakie.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "dense L2 error over the time span",
        ylabel = "solve time [s]", title = string(id, isempty(b.name) ? "" : " / " * b.name))
    for s in SOLVERS
        name = string(nameof(typeof(s[:alg])))
        # probe: one solve at the loosest tolerance with a modest step budget. Explicit
        # solvers on stiff models (lindauer, stroh) otherwise grind to maxiters at every
        # tolerance and every repeat, which took 40+ minutes and 7 GB without a result.
        kwp = b.callback === nothing ? (;) : (; callback = b.callback)
        probe = try
            solve(b.prob, s[:alg]; abstol = ABSTOLS[1], reltol = RELTOLS[1], maxiters = 10^5, kwp...)
        catch
            nothing
        end
        if probe === nothing || !successful_retcode(probe)
            push!(rows, (; solver = name, ok = false, note = "probe solve failed (stiff problem, maxiters 1e5)"))
            continue
        end
        # a solver that already needs seconds at the loosest tolerance would take hours over the
        # whole sweep (tighter tolerances x repeats); report it as too slow instead of timing it.
        # The first probe includes compilation, so time a second one.
        tprobe = @elapsed solve(b.prob, s[:alg]; abstol = ABSTOLS[1], reltol = RELTOLS[1], maxiters = 10^5, kwp...)
        if tprobe > PROBE_SECONDS
            push!(rows, (; solver = name, ok = false, note = "too slow: probe took $(round(tprobe; digits = 1)) s"))
            continue
        end
        # Walk the tolerance grid from loose to tight, one level at a time, and stop a solver
        # at the first level that fails, gives a non-finite error, or takes longer than
        # LEVEL_SECONDS. A solver that cannot keep up is reported up to where it got, instead
        # of grinding through the tight tolerances with a growing dense solution (on GitHub
        # runners that ran out of memory).
        good = Tuple{Float64, Float64}[]
        stopped = ""
        for (a, r) in zip(ABSTOLS, RELTOLS)
            local wp
            tl = @elapsed wp = try
                benchmark_problem(b, [s]; reference = refsol, abstols = [a], reltols = [r], numruns = 3, maxiters = 10^5)
            catch e
                stopped = first(sprint(showerror, e), 120)
                nothing
            end
            wp === nothing && break
            w = wp.wps[1]
            e = w.errors[1] isa NamedTuple ? w.errors[1].L2 : w.errors[1]   # DiffEqDevTools 3: NamedTuple per tolerance; dense L2
            (isfinite(e) && e > 0 && isfinite(w.times[1])) || (stopped = "non-finite error at abstol $a"; break)
            push!(good, (e, w.times[1]))
            GC.gc()
            tl > LEVEL_SECONDS && (stopped = "stopped at abstol $a: level took $(round(tl; digits = 1)) s"; break)
        end
        isempty(good) && (push!(rows, (; solver = name, ok = false, note = isempty(stopped) ? "no finite errors" : stopped)); continue)
        scatterlines!(ax, first.(good), last.(good); label = name)
        push!(rows, (; solver = name, ok = true, errors = first.(good), times = last.(good)))
    end
    axislegend(ax; position = :rt, labelsize = 9)
    save(joinpath(OUT, string(id, isempty(b.name) ? "" : "__" * replace(b.name, r"[^A-Za-z0-9_]" => "_"), ".png")), fig)
    return (; id, name = b.name, status = "ok", reference = string(nameof(typeof(refalg))), rows, scale)
end

# fastest solver whose l2 error is at most `target` (absolute)
function fastest_at(rows, target)
    best = nothing
    for r in rows
        r.ok || continue
        ts = [t for (e, t) in zip(r.errors, r.times) if e <= target]
        isempty(ts) && continue
        t = minimum(ts)
        (best === nothing || t < best[2]) && (best = (r.solver, t))
    end
    return best
end

function main(ids)
    results = []
    for id in ids
        m = Module(Symbol("Bench_", id))
        Base.include(m, joinpath(models_dir(), id, "reproduce.jl"))
        isdefined(m, :benchmark_setups) || (@warn "$id has no benchmark_setups"; continue)
        for b in Base.invokelatest(m.benchmark_setups)
            @info "benchmarking" id b.name
            push!(results, bench_one(id, b))
        end
    end
    io = IOBuffer()
    println(io, "| model | problem | reference | fastest to relative dense L2 error <= 1e-6 | time [ms] | solvers that failed |")
    println(io, "|---|---|---|---|---|---|")
    toml = Dict{String, Any}()
    for r in results
        if r.status != "ok"
            println(io, "| ", r.id, " | ", r.name, " | ", r.status, " | | | |")
            continue
        end
        best = fastest_at(r.rows, 1e-6 * r.scale)
        failed = join([x.solver for x in r.rows if !x.ok], ", ")
        println(io, "| ", r.id, " | ", r.name, " | ", r.reference, " | ", best === nothing ? "none" : best[1], " | ",
            best === nothing ? "" : @sprintf("%.3f", 1e3 * best[2]), " | ", failed, " |")
        toml[string(r.id, "/", r.name)] = Dict("reference" => r.reference, "scale" => r.scale,
            "solvers" => Dict(x.solver => (x.ok ? Dict("errors" => x.errors, "times" => x.times) : Dict("failed" => x.note)) for x in r.rows))
    end
    write(joinpath(OUT, "summary.md"), String(take!(io)))
    open(f -> TOML.print(f, toml; sorted = true), joinpath(OUT, "summary.toml"), "w")
    print(read(joinpath(OUT, "summary.md"), String))
end

main(isempty(ARGS) ? model_ids() : ARGS)
