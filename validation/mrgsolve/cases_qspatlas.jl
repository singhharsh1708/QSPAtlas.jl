# Deterministic cases 1-7 in QSPAtlas. Writes out/jl_<case>.csv (time, DEPOT, CENT, PERIPH, RESP, CP).
# Run from QSPAtlas.jl: julia --startup-file=no --threads=4 --project=. validation/mrgsolve/cases_qspatlas.jl
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
const PNAMES = Dict(:CL => sys.CL, :V1 => sys.V1, :Q => sys.Q, :V2 => sys.V2, :KA => sys.KA,
                    :KIN => sys.KIN, :KOUT => sys.KOUT, :IMAX => sys.IMAX, :IC50 => sys.IC50)
const BASE = (CL = 2.0, V1 = 10.0, Q = 3.0, V2 = 20.0, KA = 1.1, KIN = 10.0, KOUT = 0.1, IMAX = 0.8, IC50 = 2.0)

function model(p::NamedTuple)
    θ = merge(BASE, p)
    return PopulationModel(sys; routes = ROUTES, θ, ηnames = Symbol[], Ω = zeros(0, 0),
        params = (θ, η, c) -> [PNAMES[k] => θ[k] for k in keys(BASE)],
        init = (θ, η, c) -> [sys.R => θ.KIN / θ.KOUT],
        observed = (DEPOT = sys.Ad, CENT = sys.Ac, PERIPH = sys.Ap, RESP = sys.R, CP = sys.Ac / sys.V1),
        error = (DEPOT = AdditiveError(0), CENT = AdditiveError(0), PERIPH = AdditiveError(0), RESP = AdditiveError(0), CP = AdditiveError(0)),
        tspan = (0.0, 1.0), solver = Vern9(), solve_kwargs = (abstol = 1e-10, reltol = 1e-10, maxiters = 10^8))
end

const ONE = (Q = 0.0,)
ctl(; F1 = 1.0, L1 = 0.0, F2 = 1.0, L2 = 0.0) =
    DoseControl(bioav = Dict(:depot => F1, :central => F2), lag = Dict(:depot => L1, :central => L2))

cases = [
    (:c1, ONE, 96, Regimen(100; cmt = :central, ii = 12, addl = 5), ctl()),
    (:c2a, (;), 192, Regimen(100; cmt = :central, rate = 100, ii = 12, addl = 13), ctl()),
    (:c2b, (;), 192, Regimen(Regimen(100; cmt = :central, rate = 100, ii = 12, addl = 13),
                             Regimen(50; time = 0.5, cmt = :central, rate = 10, ii = 12, addl = 13)), ctl()),
    (:c2c, (;), 192, Regimen(150; cmt = :central, rate = 10, ii = 12, addl = 13), ctl()),
    (:c3a, ONE, 168, Regimen(100; cmt = :depot, ii = 24, addl = 6), ctl(F1 = 0.6, L1 = 0.75)),
    (:c3b, ONE, 110, Regimen(100; cmt = :depot, ii = 12, addl = 6), ctl(L1 = 13.0)),
    (:c4a, ONE, 72, Regimen(100; cmt = :central, rate = 20, ii = 24, addl = 2), ctl(F2 = 0.5)),
    (:c4b, ONE, 72, Regimen(100; cmt = :central, duration = 5, ii = 24, addl = 2), ctl(F2 = 0.5)),
    (:c4c, ONE, 72, Regimen(100; cmt = :central, rate = 20, ii = 24, addl = 2), ctl(F2 = 0.5, L2 = 1.5)),
    (:c5a, ONE, 48, Regimen(100; cmt = :central, ii = 12, addl = 3, ss = true), ctl()),
    (:c5b, ONE, 72, Regimen(100; cmt = :depot, ii = 24, addl = 2, ss = true), ctl(F1 = 0.6, L1 = 0.75)),
    (:c5b2, ONE, 48, Regimen(100; cmt = :depot, ii = 12, addl = 3, ss = true), ctl(L1 = 13.0)),
    (:c5c, (;), 48, Regimen(100; cmt = :central, rate = 20, ii = 12, addl = 3, ss = true), ctl()),
    (:c5d, (;), 48, Regimen(180; cmt = :central, rate = 10, ii = 12, addl = 3, ss = true), ctl()),
    (:c5e, ONE, 48, Regimen(100; cmt = :central, ii = 12, addl = 3, ss = true), ctl()),
    (:c5f, merge(ONE, (KIN = 0.05, KOUT = 0.0005)), 48, Regimen(100; cmt = :central, ii = 12, addl = 3, ss = true), ctl()),
    (:c5g, ONE, 60, Regimen(Regimen(500; cmt = :central), Regimen(100; time = 24, cmt = :central, ii = 12, addl = 2, ss = true)), ctl()),
    (:c6a, ONE, 24, Regimen(Regimen(100; cmt = :central, rate = 10), Reset(4.0), Regimen(50; time = 8, cmt = :central)), ctl()),
    (:c6b, ONE, 24, Regimen(Regimen(100; cmt = :central, rate = 10), Regimen(50; time = 4, cmt = :central, reset = true)), ctl()),
    (:c6c, ONE, 24, Regimen(Regimen(100; cmt = :central, rate = 10), Reset(4.0), Regimen(50; time = 6, cmt = :central, rate = 5)), ctl()),
    (:c6d, ONE, 24, Regimen(Regimen(100; cmt = :central), Reset(12.0)), ctl()),
    (:c6e, ONE, 12, Regimen(Regimen(100; time = 3, cmt = :depot), Reset(4.0)), ctl(L1 = 2.0)),
    (:c6f, ONE, 24, Regimen(Regimen(100; cmt = :central, ii = 6, addl = 3), Reset(8.0)), ctl()),
    (:c7a, ONE, 24, Regimen(Regimen(100; cmt = :central), Regimen(50; cmt = :central),
                            Regimen(100; time = 6, cmt = :central), Regimen(50; time = 6, cmt = :central)), ctl()),
    (:c7b, ONE, 24, Regimen(Regimen(100; cmt = :central, rate = 50), Regimen(30; time = 2, cmt = :central),
                            Regimen(60; time = 2, cmt = :central, rate = 20)), ctl()),
]

const SSO = SteadyStateOptions(rtol = 1e-10, atol = 1e-12, maxcycles = 20000)
const COLS = (:DEPOT, :CENT, :PERIPH, :RESP, :CP)
for (name, p, tend, reg, control) in cases
    grid = collect(0.0:0.25:tend)
    s = Subject("1"; regimen = reg, obstimes = grid, control)
    res = try
        zero_re_prediction(model(p), s; tspan = (0.0, Float64(tend)), steady_state = SSO)
    catch err
        println(name, " ERROR: ", sprint(showerror, err)[1:min(end, 300)])
        continue
    end
    M = hcat(grid, [[r.ipred for r in res.obs if r.endpoint === c] for c in COLS]...)
    open(joinpath(OUT, "jl_$(name).csv"), "w") do io
        println(io, "time,DEPOT,CENT,PERIPH,RESP,CP")
        writedlm(io, M, ',')
    end
    println(name, " ", size(M, 1), " rows")
end

# c5f with the default SteadyStateOptions (maxcycles 1000): expected to refuse
try
    zero_re_prediction(model(merge(ONE, (KIN = 0.05, KOUT = 0.0005))),
        Subject("1"; regimen = Regimen(100; cmt = :central, ii = 12, addl = 3, ss = true), obstimes = [0.0]); tspan = (0.0, 48.0))
    println("c5f default options: no error")
catch err
    println("c5f default options: ", sprint(showerror, err)[1:min(end, 300)])
end
