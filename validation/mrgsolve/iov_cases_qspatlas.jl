# Dose control from the model and between-occasion variability in QSPAtlas, for comparison with
# iov_cases_mrgsolve.R. One-compartment model with first-order absorption:
#   F (depot) = TVF exp(ηF) WT / 70, ALAG (depot) = TVLAG exp(ηLAG), CL = TVCL exp(κCL of the occasion)
#
#  1. Exact case (datasets/iov_exact.csv): every subject has SUPPLIED η and κ (columns ETAF, ETALAG,
#     K1, K2 = κCL of occasions 1 and 2), so individual curves can be compared record by record.
#     Occasions from the OCC column (NOCB), WT time-varying for subject 2 (the dose at 24 uses WT 95),
#     ADDL doses for subject 3, ADDL doses running into occasion 2 and WT 95 for subject 5 (they keep
#     the F and ALAG of their record), occasion labels 1, 2, 1 for subject 4 (ETALAG 0.1, so the lagged doses do not enter exactly at an observation time).
#     Writes out/jl_iov_exact.csv (ID, TIME, DEPOT, CENT, CP), one row per observation record.
#  2. Distribution case: 5,000 subjects with sampled η and κ (out/iov_pop.csv, also read by the
#     R script). Writes out/jl_iov_pop.csv (ID, TIME, CP).
# Run from QSPAtlas.jl: julia --startup-file=no --threads=4 --project=. validation/mrgsolve/iov_cases_qspatlas.jl
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq, DelimitedFiles
using ModelingToolkit: t_nounits as t, D_nounits as D

const DIR = joinpath(@__DIR__, "datasets")
const OUT = joinpath(@__DIR__, "out")
mkpath(OUT)

@parameters CL = 2.0 V = 10.0 KA = 1.1
@variables Ad(t) = 0.0 Ac(t) = 0.0
@mtkcompile sys = System([D(Ad) ~ -KA * Ad, D(Ac) ~ KA * Ad - CL / V * Ac], t)

const Θ = (CL = 2.0, V = 10.0, KA = 1.1, F = 0.6, LAG = 0.75)
const Ω = [0.09 0.0; 0.0 0.16]          # ηF, ηLAG
const ΩIOV = [0.0625;;]                 # κCL

model = PopulationModel(sys; routes = [Route(:depot, sys.Ad), Route(:central, sys.Ac)], θ = Θ,
    ηnames = [:F, :LAG], Ω, κnames = [:CL], Ω_IOV = ΩIOV,
    params = (θ, η, c, κ) -> [sys.CL => θ.CL * exp(κ.CL), sys.V => θ.V, sys.KA => θ.KA],
    dosecontrol = (θ, η, c, κ) -> (bioav = Dict(:depot => θ.F * exp(η.F) * c.WT / 70), lag = Dict(:depot => θ.LAG * exp(η.LAG))),
    observed = (DEPOT = sys.Ad, CENT = sys.Ac, CP = sys.Ac / sys.V),
    error = (DEPOT = AdditiveError(0), CENT = AdditiveError(0), CP = AdditiveError(0)),
    tspan = (0.0, 1.0), solver = Vern9(), solve_kwargs = (abstol = 1e-12, reltol = 1e-12, maxiters = 10^8))

# ---- 1. exact case -------------------------------------------------------------------------
file = joinpath(DIR, "iov_exact.csv")
data = read_nmdata(file; cmt = Dict(1 => :depot, 2 => :central), endpoints = :all, occasion = :OCC, drop = [:ETAF, :ETALAG, :K1, :K2])
raw = readdlm(file, ','; header = true)
hdr = vec(raw[2])
col(name) = raw[1][:, findfirst(==(name), hdr)]
ids = string.(Int.(col("ID")))
first_of(id, name) = col(name)[findfirst(==(id), ids)]
rows = Any[]
for s in data.subjects
    η = (F = first_of(s.id, "ETAF"), LAG = first_of(s.id, "ETALAG"))
    κ = Dict(1 => (CL = first_of(s.id, "K1"),), 2 => (CL = first_of(s.id, "K2"),))
    r = simulate(model, s; η, κ, residual = false)
    v(ep) = [o.ipred for o in r.obs if o.endpoint === ep]
    for (k, tt) in enumerate(s.obstimes)
        push!(rows, [s.id, tt, v(:DEPOT)[k], v(:CENT)[k], v(:CP)[k]])
    end
    println("subject ", s.id, ": occasions ", s.occasion_times, " => ", s.occasion_ids, ", η ", η, ", κ ", κ)
end
open(joinpath(OUT, "jl_iov_exact.csv"), "w") do io
    println(io, "ID,TIME,DEPOT,CENT,CP")
    writedlm(io, permutedims(reduce(hcat, rows)), ',')
end

# ---- 2. distribution case ----------------------------------------------------------------
const N = 5000
grid1 = [0.5, 1, 2, 4, 8, 12, 23]
popfile = joinpath(OUT, "iov_pop.csv")
open(popfile, "w") do io
    println(io, "ID,TIME,EVID,AMT,CMT,MDV,OCC,WT")
    for i in 1:N, day in 0:1
        println(io, i, ",", 24day, ",1,100,1,1,", day + 1, ",70")
        for g in grid1
            println(io, i, ",", 24day + g, ",0,.,2,0,", day + 1, ",70")
        end
    end
end
pd = read_nmdata(popfile; cmt = Dict(1 => :depot, 2 => :central), endpoints = Dict(2 => :CP), occasion = :OCC)
res = simulate(model, pd.subjects; seed = 20261105, residual = false)
open(joinpath(OUT, "jl_iov_pop.csv"), "w") do io
    println(io, "ID,TIME,CP")
    for o in res.obs
        o.endpoint === :CP && println(io, o.id, ",", o.time, ",", o.ipred)
    end
end
println("distribution case: ", length(pd.subjects), " subjects, ", count(s -> s.retcode !== :Success, res.subjects), " failed")
