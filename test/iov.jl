# Dose control from the population model (F, ALAG, modelled rate and duration with random
# effects and covariates) and between-occasion variability (IOV). Own module. Expected values are
# closed forms of a one-compartment model with an AUC state (k = CL / V):
#   bolus D with bioavailability F and lag L into the central compartment:
#     A(t) = F D exp(-k (t - L)) for t >= L,  AUC(T) = F D / CL (1 - exp(-k (T - L)))
module IOVTests

using QSPAtlas, Test
using QSPAtlas: stream_rng, sample_kappa, occasion_at, covariates_at
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using Random: Xoshiro, shuffle
using Statistics: mean, quantile, cov, cor

const TOL = (abstol = 1e-11, reltol = 1e-11)
const Z95 = 1.6448536269514722
φ(z) = exp(-z^2 / 2) / sqrt(2π)

@parameters CL = 2.0 V = 10.0 ka = 1.0 Rin = 0.0
@variables Ad(t) = 0.0 A(t) = 0.0 AUC(t) = 0.0
@mtkcompile pk = System([D(Ad) ~ -ka * Ad, D(A) ~ ka * Ad - CL / V * A + Rin, D(AUC) ~ A / V], t)
@mtkcompile pkss = System([D(Ad) ~ -ka * Ad, D(A) ~ ka * Ad - CL / V * A + Rin], t)    # no AUC: has a steady state
routes(sys) = [Route(:depot, sys.Ad), Route(:central, sys.A; rate = sys.Rin)]
const Θ = (CL = 2.0, V = 10.0, F = 0.6, LAG = 1.5, D1 = 3.0, R1 = 25.0)

# CL and V scaled by exp(κ) when the model has IOV (κ names CL, V); WT scales nothing.
params3(θ, η, c) = [CL => θ.CL, V => θ.V, ka => 1.0]
params4(θ, η, c, κ) = [CL => θ.CL * exp(get(κ, :CL, 0.0)), V => θ.V * exp(get(κ, :V, 0.0)), ka => 1.0]

mk(; ηnames = Symbol[], Ω = zeros(0, 0), κnames = Symbol[], Ω_IOV = zeros(0, 0), dosecontrol = nothing,
    params = isempty(κnames) ? params3 : params4, tspan = (0.0, 100.0), ss = false) =
    PopulationModel(ss ? pkss : pk; routes = routes(ss ? pkss : pk), θ = Θ, ηnames, Ω, κnames, Ω_IOV, params, dosecontrol,
        observed = ss ? (C = A / V, A = A) : (C = A / V, A = A, AUC = pk.AUC),
        error = ss ? (C = ProportionalError(0.1), A = ProportionalError(0.1)) :
                (C = ProportionalError(0.1), A = ProportionalError(0.1), AUC = ProportionalError(0.1)),
        tspan, solver = Tsit5(), solve_kwargs = TOL)

vals(res, ep; time = nothing, id = nothing) =
    [o.ipred for o in res.obs if o.endpoint === ep && (time === nothing || o.time == time) && (id === nothing || o.id == id)]
diag2(M) = [M[i, i] for i in 1:size(M, 1)]
relerr(a, b) = maximum(abs.(a .- b) ./ max.(abs.(b), 1e-300))

@testset "IIV on F: individual AUC exact, population AUC log-normal" begin
    ω, n = 0.3, 4000
    m = mk(; ηnames = [:F], Ω = [ω^2;;], dosecontrol = (θ, η, c) -> (bioav = Dict(:central => θ.F * exp(η.F)),))
    pop = virtual_population(n; regimen = Regimen(100.0; cmt = :central), obstimes = [100.0])
    res = simulate(m, pop; seed = 11, residual = false)
    auc = vals(res, :AUC)
    scale = 100 / Θ.CL * (1 - exp(-Θ.CL / Θ.V * 100))
    @test relerr(auc, [Θ.F * exp(s.η.F) * scale for s in res.subjects]) < 1e-8
    # log AUC ~ N(log(F scale), ω²): sample quantiles against the exact ones (order-statistic SE)
    for (q, z) in ((0.05, -Z95), (0.5, 0.0), (0.95, Z95))
        se = ω * sqrt(q * (1 - q) / n) / φ(z)
        @test abs(log(quantile(auc, q)) - (log(Θ.F * scale) + ω * z)) < 4.5se
    end
    # bioavailability on an oral dose too, and for a fixed-duration infusion it scales the rate
    m2 = mk(; ηnames = [:F], Ω = [ω^2;;], dosecontrol = (θ, η, c) -> (bioav = Dict(:depot => θ.F * exp(η.F), :central => θ.F * exp(η.F)),))
    s2 = Subject("x"; regimen = Regimen(Regimen(100.0; cmt = :depot), Regimen(50.0; time = 1.0, cmt = :central, duration = 2.0)), obstimes = [100.0])
    r2 = simulate(m2, s2; η = [0.2], residual = false)
    k = Θ.CL / Θ.V
    F = Θ.F * exp(0.2)
    @test vals(r2, :AUC)[1] ≈ F * 150 / Θ.CL rtol = 1e-7     # all of it is cleared by t = 100 (to 1e-8)
end

@testset "IIV on ALAG: individual curves exact, lag log-normal" begin
    ω, n = 0.4, 4000
    m = mk(; ηnames = [:LAG], Ω = [ω^2;;], dosecontrol = (θ, η, c) -> (lag = Dict(:central => θ.LAG * exp(η.LAG)),))
    tq = [Θ.LAG * exp(-ω * Z95), Θ.LAG, Θ.LAG * exp(ω * Z95)]           # lag quantiles 5, 50, 95 %
    pop = virtual_population(n; regimen = Regimen(100.0; cmt = :central), obstimes = vcat(tq, 6.0))
    res = simulate(m, pop; seed = 5, residual = false)
    k = Θ.CL / Θ.V
    lags = [Θ.LAG * exp(s.η.LAG) for s in res.subjects]
    a6 = vals(res, :A; time = 6.0)
    ok = lags .< 6.0
    @test relerr(a6[ok], 100 .* exp.(-k .* (6.0 .- lags[ok]))) < 1e-8
    # P(nothing given yet at t) = P(lag > t) = 0.95, 0.5, 0.05
    for (tt, p) in zip(tq, (0.95, 0.5, 0.05))
        frac = mean(vals(res, :A; time = tt) .== 0.0)
        @test abs(frac - p) < 4.5 * sqrt(p * (1 - p) / n)
    end
end

@testset "modelled rate and duration with IIV (RATE -1 / -2)" begin
    m = mk(; ηnames = [:D1, :R1], Ω = [0.09 0.0; 0.0 0.04],
        dosecontrol = (θ, η, c) -> (duration = Dict(:central => θ.D1 * exp(η.D1)), rate = Dict(:depot => θ.R1 * exp(η.R1))))
    k = Θ.CL / Θ.V
    s = Subject("d"; regimen = Regimen(100.0; cmt = :central, duration = :model), obstimes = [1.0, 20.0])
    for ηv in ([0.3, 0.0], [-0.5, 0.0])
        Di = Θ.D1 * exp(ηv[1])
        r = simulate(m, s; η = ηv, residual = false)
        inf(t) = t <= Di ? 100 / Di / k * (1 - exp(-k * t)) : 100 / Di / k * (1 - exp(-k * Di)) * exp(-k * (t - Di))
        @test vals(r, :A) ≈ inf.([1.0, 20.0]) rtol = 1e-8
    end
    # a modelled rate needs a route with a rate parameter: :depot has none, a clear error
    s2 = Subject("r"; regimen = Regimen(100.0; cmt = :depot, rate = :model), obstimes = [1.0])
    @test_throws ArgumentError simulate(m, s2; η = :zero)
end

@testset "covariate-dependent F: the dose uses the covariates at its record" begin
    m = mk(; dosecontrol = (θ, η, c) -> (bioav = Dict(:central => 0.5 * c.WT / 70),))
    k = Θ.CL / Θ.V
    reg = Regimen(Regimen(100.0; cmt = :central), Regimen(100.0; time = 24.0, cmt = :central))
    # LOCF: the change at 24 is in force for the dose record at 24
    s = Subject("l"; regimen = reg, covariates = (WT = 70.0,), timevarying = [24.0 => (WT = 140.0,)], obstimes = [12.0, 30.0])
    r = zero_re_prediction(m, s)
    @test vals(r, :A) ≈ [50exp(-12k), 50exp(-30k) + 100exp(-6k)] rtol = 1e-8
    # an ADDL dose belongs to its record: it keeps the record's F (NONMEM rule for ALAG, mrgsolve)
    sa = Subject("a"; regimen = Regimen(100.0; cmt = :central, ii = 24.0, addl = 1), covariates = (WT = 70.0,),
        timevarying = [24.0 => (WT = 140.0,)], obstimes = [30.0])
    @test vals(zero_re_prediction(m, sa), :A) ≈ [50exp(-30k) + 50exp(-6k)] rtol = 1e-8
    # NOCB dataset: WT on the dose record at 24 is 140; the observation at 12 carries WT 70
    rec(time, evid, amt, wt) = (ID = 1, TIME = time, EVID = evid, AMT = amt, CMT = 2, WT = wt)
    rows = [rec(0.0, 1, 100, 70), rec(12.0, 0, ".", 70), rec(24.0, 1, 100, 140), rec(30.0, 0, ".", 140)]
    data = read_nmdata(rows; cmt = Dict(2 => :central), endpoints = Dict(2 => :A))
    rd = simulate(m, data.subjects; η = :zero, residual = false)
    @test vals(rd, :A) ≈ [50exp(-12k), 50exp(-30k) + 100exp(-6k)] rtol = 1e-8
    # a lagged dose keeps the F of its record even if the covariate changes before it enters
    s3 = Subject("g"; regimen = Regimen(100.0; time = 20.0, cmt = :central), covariates = (WT = 70.0,),
        timevarying = [22.0 => (WT = 140.0,)], obstimes = [30.0], control = DoseControl(lag = Dict(:central => 4.0)))
    @test vals(zero_re_prediction(m, s3), :A) ≈ [50exp(-6k)] rtol = 1e-8
    # steady state with model F: trough = F D e^{-k ii} / (1 - e^{-k ii})
    F = 0.5 * 105 / 70
    mss = mk(; dosecontrol = (θ, η, c) -> (bioav = Dict(:central => 0.5 * c.WT / 70),), ss = true)
    @test vals(zero_re_prediction(mss, Subject("s"; regimen = Regimen(100.0; cmt = :central, ii = 12.0, ss = 1),
        covariates = (WT = 105.0,), obstimes = [0.0])), :A)[1] ≈ F * 100 / (1 - exp(-12k)) rtol = 1e-7
    # subject control and model dosecontrol combine; the same entry in both is an error
    sc = Subject("c"; regimen = Regimen(100.0; cmt = :central), covariates = (WT = 70.0,), obstimes = [5.0],
        control = DoseControl(lag = Dict(:central => 1.0)))
    @test vals(zero_re_prediction(m, sc), :A) ≈ [50exp(-4k)] rtol = 1e-8
    bad = Subject("b"; regimen = Regimen(100.0; cmt = :central), covariates = (WT = 70.0,), obstimes = [5.0],
        control = DoseControl(bioav = Dict(:central => 1.0)))
    @test_throws ArgumentError zero_re_prediction(m, bad)
    # malformed dosecontrol results
    @test_throws ArgumentError zero_re_prediction(mk(; dosecontrol = (θ, η, c) -> (F = 1.0,)), sc)
    @test_throws ArgumentError zero_re_prediction(mk(; dosecontrol = (θ, η, c) -> 1.0), sc)
    @test_throws ArgumentError zero_re_prediction(mk(; dosecontrol = (θ, η, c) -> (bioav = Dict(:central => -1.0),)), sc)
end

@testset "IOV: per-occasion draws (covariance, independence, streams)" begin
    Ωk = [0.04 0.012; 0.012 0.09]
    m = mk(; ηnames = [:CL], Ω = [0.09;;], κnames = [:CL, :V], Ω_IOV = Ωk)
    n = 4000
    pop = virtual_population(n; regimen = Regimen(100.0; cmt = :central), obstimes = Float64[], occasions = [0.0, 10.0, 20.0])
    res = simulate(m, pop; seed = 21, residual = false)
    K(o) = reduce(hcat, [collect(s.κ[o]) for s in res.subjects])'
    for o in 1:3
        S = cov(K(o))
        for i in 1:2, j in 1:2
            @test abs(S[i, j] - Ωk[i, j]) < 5 * sqrt((Ωk[i, i] * Ωk[j, j] + Ωk[i, j]^2) / n)
        end
        @test all(abs.(mean(K(o); dims = 1)) .< 5 .* sqrt.(diag2(Ωk) ./ n)')
    end
    # independent across occasions and of η
    for (o1, o2) in ((1, 2), (1, 3), (2, 3)), i in 1:2, j in 1:2
        @test abs(cor(K(o1)[:, i], K(o2)[:, j])) < 4.5 / sqrt(n)
    end
    ηs = [s.η.CL for s in res.subjects]
    @test all(abs(cor(ηs, K(o)[:, i])) < 4.5 / sqrt(n) for o in 1:3, i in 1:2)
    # each occasion has its own stream keyed by (seed, replicate, id, occasion)
    s1 = res.subjects[7]
    @test s1.κ[2] == sample_kappa(m, stream_rng(21, 1, s1.id, :kappa_2))
    # η is the same stream as without IOV; an extra occasion does not change the others
    m0 = mk(; ηnames = [:CL], Ω = [0.09;;])
    r0 = simulate(m0, pop[1:20]; seed = 21, residual = false)
    @test [s.η for s in r0.subjects] == [s.η for s in res.subjects[1:20]]
    @test !haskey(r0.subjects[1], :κ)
    pop4 = virtual_population(20; regimen = Regimen(100.0; cmt = :central), obstimes = Float64[], occasions = [0.0, 10.0, 20.0, 30.0])
    r4 = simulate(m, pop4; seed = 21, residual = false)
    @test all(r4.subjects[i].κ[o] == res.subjects[i].κ[o] for i in 1:20, o in 1:3)
    # a label that comes back is the same occasion
    sback = Subject("b"; regimen = Regimen(100.0; cmt = :central), obstimes = [25.0], occasions = [0.0 => 1, 10.0 => 2, 20.0 => 1])
    rb = simulate(m, [sback]; seed = 3)
    @test sort(collect(keys(rb.subjects[1].κ))) == [1, 2]
    @test occasion_at(sback, 25.0) == 1 && occasion_at(sback, 15.0) == 2 && occasion_at(sback, -1.0) == 1
    # replicates draw new κ; κ = :zero and zero_re_prediction set them to zero
    rr = simulate(m, pop[1:5]; seed = 21, replicates = 2, residual = false)
    @test rr.subjects[1].κ[1] != rr.subjects[6].κ[1]
    @test all(v -> all(iszero, v), values(simulate(m, pop[1:2]; κ = :zero).subjects[1].κ))
    @test all(v -> all(iszero, v), values(zero_re_prediction(m, pop[1]).subjects[1].κ))
    @test simulate(m, pop[1:2]; seed = 1).settings.κ === :sample
end
@testset "IOV: parameters and F switch exactly at occasion boundaries" begin
    m = mk(; κnames = [:CL, :V, :F], Ω_IOV = [0.04 0 0; 0 0.04 0; 0 0 0.04],
        dosecontrol = (θ, η, c, κ) -> (bioav = Dict(:central => θ.F * exp(κ.F)),))
    κ = Dict(1 => (CL = 0.2, V = -0.1, F = 0.1), 2 => (CL = -0.3, V = 0.25, F = -0.4))
    CLo(o) = Θ.CL * exp(κ[o].CL)
    Vo(o) = Θ.V * exp(κ[o].V)
    Fo(o) = Θ.F * exp(κ[o].F)
    k1, k2 = CLo(1) / Vo(1), CLo(2) / Vo(2)
    A10 = 100Fo(1) * exp(-10k1)
    # LOCF (start times): occasion 2 starts AT 10, an observation at 10 sees V of occasion 2
    s = Subject("o"; regimen = Regimen(100.0; cmt = :central), obstimes = [5.0, 10.0, 15.0], occasions = [0.0, 10.0])
    r = simulate(m, s; κ, residual = false)
    @test vals(r, :C) ≈ [100Fo(1) * exp(-5k1) / Vo(1), A10 / Vo(2), A10 * exp(-5k2) / Vo(2)] rtol = 1e-8
    # NOCB: the occasion listed at 10 holds over the interval ending at 10
    sn = Subject("n"; regimen = Regimen(100.0; cmt = :central), obstimes = [5.0, 10.0, 15.0], occasions = [10.0 => 1, 30.0 => 2],
        covariate_rule = :nocb)
    rn = simulate(m, sn; κ, residual = false)
    @test vals(rn, :C) ≈ [100Fo(1) * exp(-5k1) / Vo(1), A10 / Vo(1), A10 * exp(-5k2) / Vo(2)] rtol = 1e-8
    # F: the dose record at 10 is in occasion 2; a dose given at 9 with lag 2 keeps occasion 1's F
    s2 = Subject("f"; regimen = Regimen(Regimen(100.0; cmt = :central), Regimen(100.0; time = 10.0, cmt = :central)),
        obstimes = [15.0], occasions = [0.0, 10.0])
    a15 = A10 * exp(-5k2) + 100Fo(2) * exp(-5k2)
    @test vals(simulate(m, s2; κ, residual = false), :A) ≈ [a15] rtol = 1e-8
    # an ADDL dose in occasion 2 keeps the F of its record (occasion 1)
    s2a = Subject("fa"; regimen = Regimen(100.0; cmt = :central, ii = 10.0, addl = 1), obstimes = [15.0], occasions = [0.0, 10.0])
    @test vals(simulate(m, s2a; κ, residual = false), :A) ≈ [A10 * exp(-5k2) + 100Fo(1) * exp(-5k2)] rtol = 1e-8
    s3 = Subject("g"; regimen = Regimen(100.0; time = 9.0, cmt = :central), obstimes = [15.0], occasions = [0.0, 10.0],
        control = DoseControl(lag = Dict(:central => 2.0)))
    @test vals(simulate(m, s3; κ, residual = false), :A) ≈ [100Fo(1) * exp(-4k2)] rtol = 1e-8
    # a steady-state dose uses the occasion at its time
    s4 = Subject("s"; regimen = Regimen(100.0; time = 10.0, cmt = :central, ii = 12.0, ss = 1), obstimes = [10.0], occasions = [0.0, 10.0])
    mss = mk(; κnames = [:CL, :V, :F], Ω_IOV = [0.04 0 0; 0 0.04 0; 0 0 0.04], ss = true,
        dosecontrol = (θ, η, c, κ) -> (bioav = Dict(:central => θ.F * exp(κ.F)),))
    @test vals(simulate(mss, s4; κ, residual = false), :A) ≈ [100Fo(2) / (1 - exp(-12k2))] rtol = 1e-7
    # a single κ for every occasion, and a missing occasion in a Dict
    @test vals(simulate(m, s; κ = κ[1], residual = false), :C)[3] ≈ 100Fo(1) * exp(-15k1) / Vo(1) rtol = 1e-8
    @test_throws ArgumentError simulate(m, s; κ = Dict(1 => κ[1]))
    @test_throws ArgumentError simulate(mk(), s; κ = κ[1])          # model without IOV
end

@testset "IOV: validation of Ω_IOV and κnames" begin
    @test_throws ArgumentError mk(; κnames = [:CL], Ω_IOV = [0.04 0.0; 0.0 0.04])
    @test_throws ArgumentError mk(; κnames = [:CL, :V], Ω_IOV = [0.04 0.03; 0.0 0.04])
    @test_throws ArgumentError mk(; κnames = [:CL, :V], Ω_IOV = [0.04 0.1; 0.1 0.04])          # not PSD
    @test_throws ArgumentError mk(; κnames = [:CL, :V], Ω_IOV = [-0.04 0.0; 0.0 0.04])
    @test_throws ArgumentError mk(; κnames = [:CL, :CL], Ω_IOV = [0.04 0.0; 0.0 0.04])
    @test_throws ArgumentError mk(; κnames = [:CL], Ω_IOV = [NaN;;])
    err = try
        mk(; κnames = [:CL, :V], Ω_IOV = [0.04 0.1; 0.1 0.04])
    catch e
        e
    end
    @test occursin("Ω_IOV", err.msg)
    # κ declared but nothing takes it
    @test_throws ArgumentError mk(; κnames = [:CL], Ω_IOV = [0.04;;], params = params3)
    # singular Ω_IOV (perfectly correlated) is fine
    m = mk(; κnames = [:CL, :V], Ω_IOV = [0.04 0.04; 0.04 0.04])
    κ = sample_kappa(m, Xoshiro(1))
    @test κ.CL ≈ κ.V
    # occasion times must increase and labels be integers
    @test_throws ArgumentError Subject("x"; regimen = Regimen(1.0; cmt = :central), obstimes = [1.0], occasions = [5.0, 5.0])
    @test_throws ArgumentError Subject("x"; regimen = Regimen(1.0; cmt = :central), obstimes = [1.0], occasions = [0.0 => 1.5])
end

@testset "IOV: reproducible and independent of subject order" begin
    m = mk(; ηnames = [:F], Ω = [0.09;;], κnames = [:CL, :V], Ω_IOV = [0.04 0.0; 0.0 0.09],
        dosecontrol = (θ, η, c, κ) -> (bioav = Dict(:depot => θ.F * exp(η.F)), lag = Dict(:depot => 0.5 * exp(κ.CL))))
    pop = virtual_population(40; regimen = Regimen(100.0; cmt = :depot, ii = 12.0, addl = 3), obstimes = [2.0, 13.0, 30.0, 47.0],
        occasions = [0.0, 24.0])
    a = simulate(m, pop; seed = 8, replicates = 2)
    @test a == simulate(m, pop; seed = 8, replicates = 2)
    @test a != simulate(m, pop; seed = 9, replicates = 2)
    perm = shuffle(Xoshiro(2), 1:40)
    b = simulate(m, pop[perm]; seed = 8, replicates = 2, ensemble = QSPAtlas.EnsembleSerial())
    key(r) = Dict((o.replicate, o.id, o.time, o.endpoint) => (o.ipred, o.dv) for o in r.obs)
    @test key(a) == key(b)
    @test Dict((s.replicate, s.id) => s.κ for s in a.subjects) == Dict((s.replicate, s.id) => s.κ for s in b.subjects)
    # the single-subject method is reproducible from its rng
    @test simulate(m, pop[1]; rng = Xoshiro(4)) == simulate(m, pop[1]; rng = Xoshiro(4))
end

@testset "IOV: dataset with an occasion column" begin
    m = mk(; κnames = [:CL, :V, :F], Ω_IOV = [0.04 0 0; 0 0.04 0; 0 0 0.04],
        dosecontrol = (θ, η, c, κ) -> (bioav = Dict(:central => θ.F * exp(κ.F)),))
    κ = Dict(1 => (CL = 0.2, V = -0.1, F = 0.1), 2 => (CL = -0.3, V = 0.25, F = -0.4))
    rec(time, evid, amt, occ) = (ID = 1, TIME = time, EVID = evid, AMT = amt, CMT = 2, OCC = occ)
    rows = [rec(0.0, 1, 100, 1), rec(5.0, 0, ".", 1), rec(10.0, 0, ".", 1), rec(12.0, 1, 100, 2), rec(15.0, 0, ".", 2)]
    CLo(o) = Θ.CL * exp(κ[o].CL)
    Vo(o) = Θ.V * exp(κ[o].V)
    Fo(o) = Θ.F * exp(κ[o].F)
    k1, k2 = CLo(1) / Vo(1), CLo(2) / Vo(2)
    A10 = 100Fo(1) * exp(-10k1)
    # NOCB (default): occasion 2 is in force from just after the last occasion-1 record (10)
    d = read_nmdata(rows; cmt = Dict(2 => :central), endpoints = Dict(2 => :C), occasion = :OCC)
    @test isempty(d.covariates)
    s = d.subjects[1]
    @test s.occasion_times == [10.0, 15.0] && s.occasion_ids == [1, 2]
    @test occasion_at(s, 10.0) == 1 && occasion_at(s, 10.5) == 2 && occasion_at(s, 12.0) == 2
    A12 = A10 * exp(-2k2) + 100Fo(2)
    r = simulate(m, s; κ, residual = false)
    @test vals(r, :C) ≈ [100Fo(1) * exp(-5k1) / Vo(1), A10 / Vo(1), A12 * exp(-3k2) / Vo(2)] rtol = 1e-8
    # LOCF: occasion 2 starts at its first record (12)
    dl = read_nmdata(rows; cmt = Dict(2 => :central), endpoints = Dict(2 => :C), occasion = :OCC, covariate_rule = :locf)
    @test dl.subjects[1].occasion_times == [0.0, 12.0]
    A12l = A10 * exp(-2k1) + 100Fo(2)
    @test vals(simulate(m, dl.subjects[1]; κ, residual = false), :C)[3] ≈ A12l * exp(-3k2) / Vo(2) rtol = 1e-8
    # cohort simulation of a dataset samples κ per occasion label
    rc = simulate(m, d.subjects; seed = 2)
    @test sort(collect(keys(rc.subjects[1].κ))) == [1, 2]
    # without `occasion` the column stays a covariate (unchanged behaviour)
    @test read_nmdata(rows; cmt = Dict(2 => :central), endpoints = Dict(2 => :C)).covariates == [:OCC]
    # records at one time must agree; labels are integers; the column must exist
    bad = [rec(0.0, 1, 100, 1), rec(0.0, 0, ".", 2)]
    @test_throws ArgumentError read_nmdata(bad; cmt = Dict(2 => :central), occasion = :OCC)
    @test_throws ArgumentError read_nmdata([rec(0.0, 1, 100, 1.5)]; cmt = Dict(2 => :central), occasion = :OCC)
    @test_throws ArgumentError read_nmdata(rows; cmt = Dict(2 => :central), occasion = :OCCASION)
end

@testset "IOV and dose control: results do not depend on the number of threads (separate processes)" begin
    script = joinpath(mktempdir(), "run.jl")
    write(script, """
        using QSPAtlas, ModelingToolkit, OrdinaryDiffEq
        using ModelingToolkit: t_nounits as t, D_nounits as D
        @parameters CL = 1.0 V = 10.0 ka = 1.0
        @variables Ad(t) = 0.0 A(t) = 0.0
        @mtkcompile sys = System([D(Ad) ~ -ka * Ad, D(A) ~ ka * Ad - CL / V * A], t)
        m = PopulationModel(sys; routes = [Route(:depot, sys.Ad), Route(:central, sys.A)], θ = (CL = 2.0, V = 10.0, F = 0.7, LAG = 0.5),
            ηnames = [:F, :LAG], Ω = [0.09 0.0; 0.0 0.04], κnames = [:CL, :V], Ω_IOV = [0.04 0.01; 0.01 0.02],
            params = (θ, η, c, κ) -> [sys.CL => θ.CL * exp(κ.CL), sys.V => θ.V * exp(κ.V), sys.ka => 1.0],
            dosecontrol = (θ, η, c, κ) -> (bioav = Dict(:depot => θ.F * exp(η.F)), lag = Dict(:depot => θ.LAG * exp(η.LAG))),
            observed = (C = sys.A / sys.V,), error = (C = CombinedError(0.05, 0.1),),
            tspan = (0.0, 72.0), solver = Tsit5(), solve_kwargs = (abstol = 1e-10, reltol = 1e-10))
        pop = virtual_population(60; regimen = Regimen(100.0; cmt = :depot, ii = 24.0, addl = 2), obstimes = [1.0, 6.0, 25.0, 50.0, 72.0],
            occasions = [0.0, 24.0, 48.0])
        res = simulate(m, pop; seed = 99, replicates = 3)
        open(ARGS[1], "w") do io
            println(io, "threads ", Threads.nthreads())
            for o in res.obs
                println(io, o.replicate, ",", o.id, ",", o.time, ",", o.endpoint, ",", o.ipred, ",", o.dv, ",", o.retcode)
            end
            for r in res.subjects
                println(io, r.replicate, ",", r.id, ",", r.η, ",", sort(collect(r.κ); by = first), ",", r.retcode)
            end
        end
        """)
    dir = dirname(script)
    outs = Dict(n => joinpath(dir, "out$n.txt") for n in (1, 4))
    jl = Base.julia_cmd()
    proj = dirname(Base.active_project())
    procs = [run(pipeline(`$jl --startup-file=no --threads=$n --project=$proj $script $(outs[n])`; stderr = devnull); wait = false) for n in (1, 4)]
    foreach(wait, procs)
    @test all(success, procs)
    l1, l4 = readlines(outs[1]), readlines(outs[4])
    @test l1[1] == "threads 1" && l4[1] == "threads 4"
    @test length(l1) > 900 && l1[2:end] == l4[2:end]
end

@testset "IOV: pcVPC uses the typical prediction, η = 0 and κ = 0" begin
    # Bergstrand et al. 2011: PRED is the typical population prediction, so the occasion effects are
    # zero too, whatever κ the simulation itself uses
    m = mk(; ηnames = [:CL], Ω = [0.09;;], κnames = [:CL, :V], Ω_IOV = [0.04 0.0; 0.0 0.09])
    subs = [Subject("i$i"; regimen = Regimen(d; cmt = :central), obstimes = [2.0, 12.0], occasions = [0.0, 6.0])
            for (i, d) in enumerate((100.0, 100.0, 200.0, 200.0))]
    typ(s, tm) = only(o.ipred for o in zero_re_prediction(m, s).obs if o.endpoint === :C && o.time == tm)
    want = [(typ(subs[1], tm) + typ(subs[3], tm)) / 2 for tm in (2.0, 12.0)]       # median over 2 + 2 records
    for κ in (:sample, (CL = 0.5, V = -0.3))
        v = vpc(m, subs; endpoint = :C, q = (0.5,), replicates = 4, seed = 3, pc = true, κ)
        @test v.pred_bin ≈ want rtol = 1e-9
    end
    # the SimResult method: predictions with random occasion effects are refused, typical ones match the model method
    res = simulate(m, subs; replicates = 4, seed = 3, κ = :sample)
    @test_throws ArgumentError vpc(res; endpoint = :C, q = (0.5,), pc = true, pred = simulate(m, subs; η = :zero, κ = :sample, residual = false))
    @test_throws ArgumentError vpc(res; endpoint = :C, q = (0.5,), pc = true, pred = simulate(m, subs; residual = false))
    typical = simulate(m, subs; η = :zero, κ = :zero, residual = false)
    @test isequal(vpc(res; endpoint = :C, q = (0.5,), pc = true, pred = typical).rows,
        vpc(m, subs; endpoint = :C, q = (0.5,), replicates = 4, seed = 3, pc = true, κ = :sample).rows)
end

end # module IOVTests
