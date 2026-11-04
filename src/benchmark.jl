using DiffEqDevTools: DiffEqDevTools, WorkPrecisionSet, TestSolution
using SciMLBase: solve

"""
    BenchmarkSetup(prob; name = "", callback = nothing, saveat = nothing)

What a model exposes for solver benchmarking: the `ODEProblem` used in its
reproduction, any dosing `callback`, and the `saveat` grid the figures need.
Each model's `protocol.jl` defines `benchmark_setups()` returning a vector of these.
"""
Base.@kwdef struct BenchmarkSetup
    name::String = ""
    prob::Any
    callback::Any = nothing
    saveat::Any = nothing
end
BenchmarkSetup(prob; kw...) = BenchmarkSetup(; prob, kw...)

"""
    benchmark_problem(setup, solvers; reference_alg = nothing, reference = nothing, abstols, reltols, numruns = 10, kwargs...)

Work-precision set for one model problem. Solves a reference with `reference_alg` at
abstol = reltol = 1e-12 (or uses `reference`, a solution computed that way beforehand), then times every solver in `solvers` (DiffEqDevTools
format, e.g. `[Dict(:alg => Tsit5())]`) across the tolerance grid.

Error is the dense L2 error (`error_estimate = :L2`, `dense_errors = true`): both
solutions are compared through their interpolants over the whole time span. Sparse
errors at saved points are wrong for dosed models, because a dosing callback saves
each dose time twice (before and after the jump) and the comparison then measures the
dose size instead of the solver error. `saveat` is therefore not passed to the timed
solves; timings include dense output.
"""
function benchmark_problem(b::BenchmarkSetup, solvers; reference_alg = nothing, reference = nothing,
        abstols = 1.0 ./ 10.0 .^ (5:10), reltols = 1.0 ./ 10.0 .^ (2:7), numruns = 10, kwargs...)
    kw = b.callback === nothing ? (;) : (; callback = b.callback)
    reference === nothing && reference_alg === nothing && throw(ArgumentError("give reference_alg or a precomputed reference solution"))
    ref = reference === nothing ? solve(b.prob, reference_alg; abstol = 1e-12, reltol = 1e-12, kw...) : reference
    return WorkPrecisionSet(b.prob, abstols, reltols, solvers; appxsol = TestSolution(ref), error_estimate = :L2,
        dense_errors = true, numruns, names = [string(nameof(typeof(s[:alg]))) for s in solvers], kw..., kwargs...)
end
