# Case 8, sampling check without the ODE solver: 10^6 subjects' eta and residual draws taken
# exactly as `simulate` takes them (stream_rng(seed, replicate, id, :eta / :residual),
# sample_eta, sample_dv), pushed through the closed-form one-compartment oral solution.
# Writes out/rng_check_jl.csv (time, var, q, value) for pop_rng_check.R.
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq, Statistics, LinearAlgebra, DelimitedFiles
using QSPAtlas: stream_rng, sample_eta, sample_dv
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters CL = 5.0 V = 50.0
@variables A(t) = 0.0
@mtkcompile sys = System([D(A) ~ -CL / V * A], t)
Ω = [0.09 0.024; 0.024 0.04]
m = PopulationModel(sys; routes = [Route(:central, sys.A)], θ = (CL = 5.0, V = 50.0), ηnames = [:CL, :V], Ω,
    params = (θ, η, c) -> [sys.CL => θ.CL * exp(η.CL), sys.V => θ.V * exp(η.V)],
    observed = (CP = sys.A / sys.V,), error = (CP = ProportionalError(0.15),), tspan = (0.0, 1.0), solver = Tsit5())

const N = 10^6
const OBST = [0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 24.5, 26, 30, 36, 48, 50, 54, 60, 72]
ipred(t, CLi, Vi) = (k = CLi / Vi; sum(t >= td ? 100 / (Vi * (1 - k)) * (exp(-k * (t - td)) - exp(-(t - td))) : 0.0 for td in (0.0, 24.0, 48.0)))
E = zeros(N, 2)
IP = zeros(N, length(OBST)); DV = similar(IP)
err = ProportionalError(0.15)
Threads.@threads for i in 1:N
    η = sample_eta(m, stream_rng(11, 1, string(i), :eta))
    E[i, 1] = η.CL; E[i, 2] = η.V
    r = stream_rng(11, 1, string(i), :residual)
    for (j, tt) in enumerate(OBST)
        IP[i, j] = ipred(tt, 5 * exp(η.CL), 50 * exp(η.V))
        DV[i, j] = sample_dv(err, IP[i, j], r)
    end
end
println("eta mean ", round.(vec(mean(E; dims = 1)); sigdigits = 3), " (SE ~ ", round(0.3 / sqrt(N); sigdigits = 2), ")")
println("eta cov  ", round.(cov(E); digits = 5), " target ", Ω)
rows = Any[]
for (j, tt) in enumerate(OBST), (v, X) in (("IPRED", IP), ("DV", DV)), q in (0.05, 0.5, 0.95)
    push!(rows, (tt, v, q, quantile(view(X, :, j), q)))
end
open(joinpath(@__DIR__, "out", "rng_check_jl.csv"), "w") do io
    println(io, "time,var,q,value")
    writedlm(io, rows, ',')
end
