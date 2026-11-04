# Case 8 in QSPAtlas: one-compartment oral, log-normal IIV on CL and V (SD 0.3 / 0.2, corr 0.4),
# proportional error 0.15, 5000 subjects. Writes out/pop_jl.csv, out/pop_jl_eta.csv,
# out/pop_jl_zero.csv (zero-RE prediction and the same subject in a model without random effects).
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq, DelimitedFiles
using ModelingToolkit: t_nounits as t, D_nounits as D

const OUT = joinpath(@__DIR__, "out")
@parameters CL = 5.0 V = 50.0 KA = 1.0
@variables Ad(t) = 0.0 Ac(t) = 0.0
@mtkcompile sys = System([D(Ad) ~ -KA * Ad, D(Ac) ~ KA * Ad - CL / V * Ac], t)

const OBST = [0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 24.5, 26, 30, 36, 48, 50, 54, 60, 72]
const REG = Regimen(100; cmt = :depot, ii = 24, addl = 2)
const KW = (abstol = 1e-10, reltol = 1e-10)
const Ω = [0.09 0.4*0.3*0.2; 0.4*0.3*0.2 0.04]

popmodel(ηnames, Ω) = PopulationModel(sys; routes = [Route(:depot, sys.Ad)], θ = (CL = 5.0, V = 50.0, KA = 1.0),
    ηnames, Ω,
    params = (θ, η, c) -> [sys.CL => θ.CL * exp(get(η, :CL, 0.0)), sys.V => θ.V * exp(get(η, :V, 0.0)), sys.KA => θ.KA],
    observed = (CP = sys.Ac / sys.V,), error = (CP = ProportionalError(0.15),),
    tspan = (0.0, 72.0), solver = Vern9(), solve_kwargs = KW)

m = popmodel([:CL, :V], Ω)
pop = virtual_population(5000; regimen = REG, obstimes = OBST)
res = simulate(m, pop; seed = 20261104)
open(joinpath(OUT, "pop_jl.csv"), "w") do io
    println(io, "id,time,IPRED,DV")
    writedlm(io, [(r.id, r.time, r.ipred, r.dv) for r in res.obs], ',')
end
open(joinpath(OUT, "pop_jl_eta.csv"), "w") do io
    println(io, "id,ETA1,ETA2")
    writedlm(io, [(s.id, s.η.CL, s.η.V) for s in res.subjects], ',')
end

# a second, independent run: 4 replicates of the cohort (20000 subjects) with another seed
res4 = simulate(m, pop; seed = 7, replicates = 4)
open(joinpath(OUT, "pop_jl_20k.csv"), "w") do io
    println(io, "id,time,IPRED,DV")
    writedlm(io, [(string(r.replicate, "_", r.id), r.time, r.ipred, r.dv) for r in res4.obs], ',')
end

# zero-RE prediction vs the same subject in a model with no random effects at all
grid = collect(0.0:0.25:72.0)
s = Subject("z"; regimen = REG, obstimes = grid)
z = [r.ipred for r in zero_re_prediction(m, s).obs]
det = [r.ipred for r in zero_re_prediction(popmodel(Symbol[], zeros(0, 0)), s).obs]
open(joinpath(OUT, "pop_jl_zero.csv"), "w") do io
    println(io, "time,ZERO,DET")
    writedlm(io, hcat(grid, z, det), ',')
end
println("subjects ", length(res.subjects), ", rows ", length(res.obs), ", zero-RE vs deterministic max |diff| = ", maximum(abs.(z .- det)))
