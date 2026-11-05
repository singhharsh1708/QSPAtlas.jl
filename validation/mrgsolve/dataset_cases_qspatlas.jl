# Dataset cases d1 to d7 in QSPAtlas: the NONMEM-style files in datasets/ are read with
# read_nmdata and simulated without random effects. Writes out/jl_<case>.csv with one row per
# observation record, in record order (ID, TIME, DEPOT, CENT, PERIPH, RESP, CP).
# Run from QSPAtlas.jl: julia --startup-file=no --threads=4 --project=. validation/mrgsolve/dataset_cases_qspatlas.jl
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq, DelimitedFiles
using ModelingToolkit: t_nounits as t, D_nounits as D

const DIR = joinpath(@__DIR__, "datasets")
const OUT = joinpath(@__DIR__, "out")
mkpath(OUT)

# The model of cases_qspatlas.jl, with clearance and central volume scaled by the covariate WT.
@parameters CL = 2.0 V1 = 10.0 Q = 3.0 V2 = 20.0 KA = 1.1 KIN = 10.0 KOUT = 0.1 IMAX = 0.8 IC50 = 2.0 Rd = 0.0 Rc = 0.0
@variables Ad(t) = 0.0 Ac(t) = 0.0 Ap(t) = 0.0 R(t)
eqs = [D(Ad) ~ -KA * Ad + Rd,
       D(Ac) ~ KA * Ad - (CL + Q) / V1 * Ac + Q / V2 * Ap + Rc,
       D(Ap) ~ Q / V1 * Ac - Q / V2 * Ap,
       D(R) ~ KIN * (1 - IMAX * (Ac / V1) / (IC50 + Ac / V1)) - KOUT * R]
@mtkcompile sys = System(eqs, t)
const ROUTES = [Route(:depot, sys.Ad; rate = sys.Rd), Route(:central, sys.Ac; rate = sys.Rc)]
const BASE = (CL = 2.0, V1 = 10.0, Q = 3.0, V2 = 20.0, KA = 1.1, KIN = 10.0, KOUT = 0.1, IMAX = 0.8, IC50 = 2.0)
const ENDPOINTS = (:DEPOT, :CENT, :PERIPH, :RESP, :CP)

function model(p::NamedTuple)
    θ = merge(BASE, p)
    return PopulationModel(sys; routes = ROUTES, θ, ηnames = Symbol[], Ω = zeros(0, 0),
        params = (θ, η, c) -> [sys.CL => θ.CL * (c.WT / 70)^0.75, sys.V1 => θ.V1 * c.WT / 70, sys.Q => θ.Q, sys.V2 => θ.V2,
                               sys.KA => θ.KA, sys.KIN => θ.KIN, sys.KOUT => θ.KOUT, sys.IMAX => θ.IMAX, sys.IC50 => θ.IC50],
        init = (θ, η, c) -> [sys.R => θ.KIN / θ.KOUT],
        observed = (DEPOT = sys.Ad, CENT = sys.Ac, PERIPH = sys.Ap, RESP = sys.R, CP = sys.Ac / sys.V1),
        error = NamedTuple{ENDPOINTS}(Tuple(AdditiveError(0) for _ in ENDPOINTS)),
        tspan = (0.0, 1.0), solver = Vern9(), solve_kwargs = (abstol = 1e-10, reltol = 1e-10, maxiters = 10^8))
end

const SSO = SteadyStateOptions(rtol = 1e-10, atol = 1e-12, maxcycles = 20000)
const CTL6 = DoseControl(bioav = Dict(:central => 0.5), rate = Dict(:central => 20.0), duration = Dict(:central => 5.0))
const CTL5 = DoseControl(rate = Dict(:central => 30.0))

# (case, file, parameters, covariate rule, dose control)
cases = [
    (:d1_order, "d1_order.csv", (;), :nocb, DoseControl()),
    (:d2_cov_nocb, "d2_cov.csv", (;), :nocb, DoseControl()),
    (:d2_cov_locf, "d2_cov.csv", (;), :locf, DoseControl()),
    (:d3_evid2, "d3_evid2.csv", (;), :nocb, DoseControl()),
    (:d4_ss2, "d4_ss2.csv", (;), :nocb, DoseControl()),
    (:d4_ss2_lin, "d4_ss2.csv", (IMAX = 0.0,), :nocb, DoseControl()),
    (:d5_ssinf, "d5_ssinf.csv", (;), :nocb, CTL5),
    (:d6_rate, "d6_rate.csv", (;), :nocb, CTL6),
    (:d7_reset, "d7_reset.csv", (;), :nocb, DoseControl()),
]

for (name, file, p, rule, control) in cases
    data = read_nmdata(joinpath(DIR, file); cmt = Dict(1 => :depot, 2 => :central), endpoints = :all,
        covariate_rule = rule, control)
    res = simulate(model(p), data.subjects; η = :zero, residual = false, steady_state = SSO)
    # rows of one subject and endpoint are in the subject's observation order = record order
    vals = Dict{Tuple{String, Symbol}, Vector{Float64}}()
    for o in res.obs
        push!(get!(vals, (o.id, o.endpoint), Float64[]), o.ipred)
    end
    seen = Dict{String, Int}()
    M = Matrix{Any}(undef, length(data.obs), 2 + length(ENDPOINTS))
    for (i, o) in enumerate(data.obs)
        k = seen[o.id] = get(seen, o.id, 0) + 1
        M[i, 1] = o.id
        M[i, 2] = o.time
        for (e, ep) in enumerate(ENDPOINTS)
            M[i, 2 + e] = vals[(o.id, ep)][k]
        end
    end
    open(joinpath(OUT, "jl_$(name).csv"), "w") do io
        println(io, "ID,TIME,", join(ENDPOINTS, ","))
        writedlm(io, M, ',')
    end
    println(name, ": ", size(M, 1), " observation records, order ", join([o.order === :pre ? "pre" : "post" for o in data.obs], " "))
end
