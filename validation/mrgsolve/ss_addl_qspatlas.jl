# Case d8_ss_addl in QSPAtlas: an SS=1 record after a record with ADDL doses still to come
# (datasets/d8_ss_addl.csv; subjects described in ss_addl_mrgsolve.R). The SS=1 record cancels
# the pending ADDL doses of the earlier record (NONMEM help, SS). Also writes the same subjects
# with the cancellation written out as explicit dose records (out/jl_ss_addl_explicit.csv),
# which must agree with the dataset result, and with those doses kept as dose records of their
# own (out/jl_ss_addl_kept.csv). Model, settings and output format as in
# dataset_cases_qspatlas.jl. Writes out/jl_d8_ss_addl.csv; compare with dataset_compare.R.
# Run from QSPAtlas.jl: julia --startup-file=no --threads=4 --project=. validation/mrgsolve/ss_addl_qspatlas.jl
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq, DelimitedFiles
using ModelingToolkit: t_nounits as t, D_nounits as D

const OUT = joinpath(@__DIR__, "out")
mkpath(OUT)

@parameters CL = 2.0 V1 = 10.0 Q = 3.0 V2 = 20.0 KA = 1.1 KIN = 10.0 KOUT = 0.1 IMAX = 0.8 IC50 = 2.0 Rd = 0.0 Rc = 0.0
@variables Ad(t) = 0.0 Ac(t) = 0.0 Ap(t) = 0.0 R(t)
eqs = [D(Ad) ~ -KA * Ad + Rd,
       D(Ac) ~ KA * Ad - (CL + Q) / V1 * Ac + Q / V2 * Ap + Rc,
       D(Ap) ~ Q / V1 * Ac - Q / V2 * Ap,
       D(R) ~ KIN * (1 - IMAX * (Ac / V1) / (IC50 + Ac / V1)) - KOUT * R]
@mtkcompile sys = System(eqs, t)
const ROUTES = [Route(:depot, sys.Ad; rate = sys.Rd), Route(:central, sys.Ac; rate = sys.Rc)]
const ENDPOINTS = (:DEPOT, :CENT, :PERIPH, :RESP, :CP)
const M = PopulationModel(sys; routes = ROUTES, θ = (;), ηnames = Symbol[], Ω = zeros(0, 0),
    params = (θ, η, c) -> [sys.CL => 2.0 * (c.WT / 70)^0.75, sys.V1 => 10.0 * c.WT / 70],
    init = (θ, η, c) -> [sys.R => 100.0],
    observed = (DEPOT = sys.Ad, CENT = sys.Ac, PERIPH = sys.Ap, RESP = sys.R, CP = sys.Ac / sys.V1),
    error = NamedTuple{ENDPOINTS}(Tuple(AdditiveError(0) for _ in ENDPOINTS)),
    tspan = (0.0, 1.0), solver = Vern9(), solve_kwargs = (abstol = 1e-10, reltol = 1e-10, maxiters = 10^8))
const SSO = SteadyStateOptions(rtol = 1e-10, atol = 1e-12, maxcycles = 20000)

function write_rows(path, ids, res)
    vals = Dict{Tuple{String, Symbol}, Vector{Tuple{Float64, Float64}}}()
    for o in res.obs
        push!(get!(vals, (o.id, o.endpoint), Tuple{Float64, Float64}[]), (o.time, o.ipred))
    end
    open(path, "w") do io
        println(io, "ID,TIME,", join(ENDPOINTS, ","))
        for id in ids, k in eachindex(vals[(id, :CENT)])
            tt = vals[(id, :CENT)][k][1]
            println(io, id, ",", tt, ",", join((vals[(id, ep)][k][2] for ep in ENDPOINTS), ","))
        end
    end
end

data = read_nmdata(joinpath(@__DIR__, "datasets", "d8_ss_addl.csv"); cmt = Dict(1 => :depot, 2 => :central), endpoints = :all)
res = simulate(M, data.subjects; η = :zero, residual = false, steady_state = SSO)
write_rows(joinpath(OUT, "jl_d8_ss_addl.csv"), ["1", "2", "3"], res)

# the same with the cancelled ADDL doses left out by hand
ot(s) = unique(s.obstimes)
cov = (WT = 70.0,)
s1 = Subject("1"; regimen = Regimen(Regimen(100.0; cmt = :central, ii = 12.0, addl = 2),
        Regimen(50.0; time = 30.0, cmt = :depot, ii = 24.0, addl = 1, ss = true)), obstimes = ot(data.subjects[1]), covariates = cov)
s3 = Subject("3"; regimen = Regimen(Regimen(100.0; cmt = :central, rate = 10.0, ii = 12.0, addl = 1),
        Regimen(40.0; time = 18.0, cmt = :central, ii = 12.0, ss = true)), obstimes = ot(data.subjects[3]), covariates = cov)
res2 = simulate(M, [s1, data.subjects[2], s3]; η = :zero, residual = false, steady_state = SSO)
write_rows(joinpath(OUT, "jl_ss_addl_explicit.csv"), ["1", "2", "3"], res2)
# and with the ADDL doses after the SS=1 time written as dose records of their own (record time =
# dose time, so not cancelled): what mrgsolve 2.0.1 computes (compared in ss_addl_mrgsolve.R)
k1 = Subject("1"; regimen = Regimen(Regimen(100.0; cmt = :central, ii = 12.0, addl = 2),
        Regimen(50.0; time = 30.0, cmt = :depot, ii = 24.0, addl = 1, ss = true),
        (Regimen(100.0; time = x, cmt = :central) for x in (36.0, 48.0, 60.0))...), obstimes = ot(data.subjects[1]), covariates = cov)
k3 = Subject("3"; regimen = Regimen(Regimen(100.0; cmt = :central, rate = 10.0, ii = 12.0, addl = 1),
        Regimen(40.0; time = 18.0, cmt = :central, ii = 12.0, ss = true),
        (Regimen(100.0; time = x, cmt = :central, rate = 10.0) for x in (24.0, 36.0))...), obstimes = ot(data.subjects[3]), covariates = cov)
res3 = simulate(M, [k1, data.subjects[2], k3]; η = :zero, residual = false, steady_state = SSO)
write_rows(joinpath(OUT, "jl_ss_addl_kept.csv"), ["1", "2", "3"], res3)

a = Float64.(readdlm(joinpath(OUT, "jl_d8_ss_addl.csv"), ',', skipstart = 1)[:, 3:end])
b = Float64.(readdlm(joinpath(OUT, "jl_ss_addl_explicit.csv"), ',', skipstart = 1)[:, 3:end])
println("dataset vs explicit doses: largest relative difference ", maximum(abs.(a .- b) ./ max.(abs.(b), 1e-12)))
