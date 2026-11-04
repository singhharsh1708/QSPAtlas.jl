# Tests for the population layer (subjects, population model, error models, simulate).
# The structural model is a one-compartment IV bolus, A' = -CL/V A, C = A/V, whose
# population behaviour has closed forms.
module PopulationTests

using QSPAtlas, Test
using ModelingToolkit, OrdinaryDiffEq, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D
using Random: Xoshiro, shuffle
using Statistics: mean, std, median, quantile, cov
using LinearAlgebra: eigvals, Symmetric

@parameters CL = 1.0 V = 10.0
@variables A(t) = 0.0
@mtkcompile onecmt = System([D(A) ~ -CL / V * A], t)

const TOL = (abstol = 1e-10, reltol = 1e-10)
const Z95 = 1.6448536269514722        # normal 95th percentile
φ(z) = exp(-z^2 / 2) / sqrt(2π)

# CL = θ exp(η) (WT/70)^0.75, V fixed. Amount A is dosed in the :central route.
function onecmt_model(; θCL = 2.0, V0 = 10.0, Ω = [0.09;;], ηnames = [:CL], error = ProportionalError(0.1),
        solve_kwargs = TOL, tspan = (0.0, 24.0))
    sys = onecmt
    return PopulationModel(sys;
        routes = [Route(:central, sys.A)],
        θ = (CL = θCL, V = V0),
        ηnames, Ω,
        params = (θ, η, c) -> [sys.CL => θ.CL * exp(get(η, :CL, 0.0)) * (get(c, :WT, 70.0) / 70)^0.75, sys.V => θ.V],
        observed = (C = sys.A / sys.V,),
        error = (C = error,),
        tspan, solver = Tsit5(), solve_kwargs)
end

dose100(; kw...) = Regimen(100.0; cmt = :central, kw...)
rowsof(res; id, rep = 1, endpoint = :C) = [r for r in res.obs if r.id == id && r.replicate == rep && r.endpoint == endpoint]

@testset "eta sampling" begin
    Ω = [0.09 0.03 0.0; 0.03 0.04 0.0; 0.0 0.0 0.0]       # third eta has zero variance
    m = onecmt_model(; Ω, ηnames = [:CL, :V, :Z])
    L = QSPAtlas.eta_factor(Ω)
    @test L * L' ≈ Ω atol = 1e-14
    rng = Xoshiro(2024)
    n = 40_000
    draws = [collect(QSPAtlas.sample_eta(m, rng)) for _ in 1:n]
    X = reduce(hcat, draws)'
    @test all(X[:, 3] .== 0.0)                               # exactly zero, not just small
    for i in 1:2
        @test abs(mean(X[:, i])) < 4 * sqrt(Ω[i, i] / n)
    end
    S = cov(X)
    for i in 1:2, j in 1:2
        se = sqrt((Ω[i, i] * Ω[j, j] + Ω[i, j]^2) / n)
        @test abs(S[i, j] - Ω[i, j]) < 5 * se
    end
    # draws are NamedTuples with the right names; a singular (rank one) PSD matrix is fine
    @test keys(QSPAtlas.sample_eta(m, rng)) == (:CL, :V, :Z)
    ms = onecmt_model(; Ω = [1.0 1.0; 1.0 1.0], ηnames = [:CL, :V])
    e = QSPAtlas.sample_eta(ms, rng)
    @test e.CL ≈ e.V atol = 1e-12
    # a genuinely singular PSD matrix (rank 2 of 3) keeps its variances and covariances
    B = [1.0 0.0; 0.5 1.0; 2.0 -1.0]
    Ωs = 1e-3 .* (B * B')
    @test minimum(eigvals(Symmetric(Ωs))) < 1e-15
    msing = onecmt_model(; Ω = Ωs, ηnames = [:CL, :V, :Z])
    @test QSPAtlas.eta_factor(Ωs) * QSPAtlas.eta_factor(Ωs)' ≈ Ωs rtol = 1e-8
    Xs = reduce(hcat, [collect(QSPAtlas.sample_eta(msing, rng)) for _ in 1:n])'
    Ss = cov(Xs)
    for i in 1:3, j in 1:3
        @test abs(Ss[i, j] - Ωs[i, j]) < 5 * sqrt((Ωs[i, i] * Ωs[j, j] + Ωs[i, j]^2) / n)
    end
    # the draw count does not depend on the structure of Ω
    r1, r2 = Xoshiro(5), Xoshiro(5)
    QSPAtlas.sample_eta(m, r1)
    QSPAtlas.sample_eta(onecmt_model(; Ω = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], ηnames = [:CL, :V, :Z]), r2)
    @test rand(r1) == rand(r2)
end

@testset "zero-random-effect prediction equals a direct solve" begin
    m = onecmt_model()
    obstimes = [0.5, 3.0, 7.0, 13.5, 24.0]
    s = Subject("a"; regimen = dose100(ii = 6.0, addl = 3), obstimes)
    r = zero_re_prediction(m, s)
    @test length(r.obs) == 5 && all(o -> o.retcode === :Success, r.obs)
    @test all(o -> o.dv == o.ipred, r.obs)                   # no residual error
    @test r.subjects[1].η == (CL = 0.0,)
    # direct deterministic solve with the existing bolus callback
    prob = ODEProblem(onecmt, [CL => 2.0, V => 10.0], (0.0, 24.0))
    cb = bolus_callback(onecmt, onecmt.A, repeated_doses(0.0, 6.0, 4), 100.0)
    sol = solve(prob, Tsit5(); callback = cb, abstol = 1e-10, reltol = 1e-10)
    direct = [sol(tt; idxs = onecmt.A) / 10 for tt in obstimes]      # none of these is a dose time
    @test [o.ipred for o in r.obs] ≈ direct rtol = 1e-6
    # and the closed form
    exact(tt) = sum(10.0 * exp(-0.2 * (tt - td)) for td in 0:6:18 if td <= tt)
    @test [o.ipred for o in r.obs] ≈ exact.(obstimes) rtol = 1e-7
    # observations at a dose time see the state after the dose
    s2 = Subject("b"; regimen = dose100(ii = 6.0, addl = 1), obstimes = [0.0, 6.0, 6.5])
    r2 = zero_re_prediction(m, s2)
    @test [o.ipred for o in r2.obs] ≈ [10.0, 10.0 + 10.0 * exp(-1.2), 10.0 * exp(-0.2 * 0.5) + 10.0 * exp(-0.2 * 6.5)] rtol = 1e-7
    # a fixed eta is used as given
    r3 = simulate(m, s2; η = (CL = log(2.0),), residual = false)
    @test r3.obs[1].ipred ≈ 10.0
    @test r3.obs[3].ipred ≈ 10.0 * exp(-0.4 * 0.5) + 10.0 * exp(-0.4 * 6.5) rtol = 1e-7
    @test simulate(m, s2; η = [log(2.0)], residual = false).obs == r3.obs
end

@testset "population percentiles match the log-normal closed form" begin
    θCL, V0, ω, Dose = 2.0, 10.0, 0.3, 100.0
    m = onecmt_model(; θCL, V0, Ω = [ω^2;;])
    n = 4000
    times = [2.0, 5.0, 10.0]
    pop = virtual_population(n; regimen = dose100(), obstimes = times)
    res = simulate(m, pop; seed = 11, residual = false)
    @test length(res.obs) == 3n
    for tt in times, p in (0.05, 0.5, 0.95)
        c = [o.ipred for o in res.obs if o.time == tt]
        sim_q = quantile(c, p)
        # C decreases with CL, so the p-quantile of C uses the (1-p)-quantile of CL
        z = p == 0.5 ? 0.0 : (p == 0.95 ? -Z95 : Z95)
        CLq = θCL * exp(ω * z)
        exact_q = Dose / V0 * exp(-CLq * tt / V0)
        # standard error of the sample quantile, carried through the monotone map
        se_eta = ω * sqrt(p * (1 - p) / n) / φ(z)
        tol = 4 * se_eta * CLq * tt / V0
        @test abs(log(sim_q) - log(exact_q)) < tol
    end
end

@testset "error models" begin
    n = 300_000
    ip = 5.0
    function draws(err, ipred = ip; seed = 1)
        rng = Xoshiro(seed)
        return [QSPAtlas.sample_dv(err, ipred, rng) for _ in 1:n]
    end
    # additive: dv = ipred + σ ε
    x = draws(AdditiveError(0.7))
    @test abs(mean(x) - ip) < 4 * 0.7 / sqrt(n)
    @test std(x) ≈ 0.7 rtol = 0.01
    @test QSPAtlas.sd(AdditiveError(0.7), ip) == 0.7
    # proportional: dv = ipred (1 + σ ε)
    x = draws(ProportionalError(0.2))
    @test abs(mean(x) - ip) < 4 * 0.2 * ip / sqrt(n)
    @test std(x) ≈ 0.2 * ip rtol = 0.01
    @test QSPAtlas.sd(ProportionalError(0.2), ip) == 0.2 * ip
    # combined: sd = sqrt(σadd² + (σprop ipred)²)
    c = CombinedError(0.5, 0.1)
    x = draws(c)
    want = sqrt(0.5^2 + (0.1 * ip)^2)
    @test QSPAtlas.sd(c, ip) ≈ want
    @test std(x) ≈ want rtol = 0.01
    @test abs(mean(x) - ip) < 4 * want / sqrt(n)
    # exponential: median-centred, log dv ~ N(log ipred, σ²)
    e = ExponentialError(0.3)
    x = draws(e)
    @test median(x) ≈ ip rtol = 0.01
    @test mean(log.(x)) ≈ log(ip) atol = 4 * 0.3 / sqrt(n)
    @test std(log.(x)) ≈ 0.3 rtol = 0.01
    @test mean(x) ≈ ip * exp(0.3^2 / 2) rtol = 0.01
    @test std(x) ≈ QSPAtlas.sd(e, ip) rtol = 0.02
    @test all(>(0), x)
    # negative observations are never clipped
    x = draws(AdditiveError(1.0), 0.1)
    @test any(<(0), x)
    # every draw consumes exactly one normal, even with σ = 0 or a NaN prediction
    for err in (AdditiveError(0.0), ProportionalError(0.0), CombinedError(0.0, 0.0), ExponentialError(0.0))
        r1, r2 = Xoshiro(3), Xoshiro(3)
        QSPAtlas.sample_dv(err, NaN, r1)
        randn(r2)
        @test rand(r1) == rand(r2)
    end
    @test_throws ArgumentError AdditiveError(-1.0)
    @test_throws ArgumentError ProportionalError(NaN)
    @test_throws ArgumentError CombinedError(0.1, -0.1)

    # end to end: proportional error through simulate()
    m = onecmt_model(; error = ProportionalError(0.25))
    pop = virtual_population(4000; regimen = dose100(), obstimes = [4.0])
    res = simulate(m, pop; seed = 5, η = :zero)
    rel = [(o.dv - o.ipred) / o.ipred for o in res.obs]
    @test all(o -> o.ipred == res.obs[1].ipred, res.obs)
    @test abs(mean(rel)) < 4 * 0.25 / sqrt(4000)
    @test std(rel) ≈ 0.25 rtol = 0.04
end

@testset "reproducibility" begin
    m = PopulationModel(onecmt; routes = [Route(:central, onecmt.A)], θ = (CL = 2.0, V = 10.0), ηnames = [:CL, :V],
        Ω = [0.09 0.0; 0.0 0.04],
        params = (θ, η, c) -> [onecmt.CL => θ.CL * exp(η.CL), onecmt.V => θ.V * exp(η.V)],
        observed = (C = onecmt.A / onecmt.V,), error = (C = CombinedError(0.05, 0.1),),
        tspan = (0.0, 24.0), solver = Tsit5(), solve_kwargs = TOL)
    pop = virtual_population(40; regimen = dose100(), obstimes = [1.0, 4.0, 12.0])
    a = simulate(m, pop; seed = 7, replicates = 2)
    @test a == simulate(m, pop; seed = 7, replicates = 2)                      # same seed, same result
    @test a == simulate(m, pop; seed = 7, replicates = 2, ensemble = SciMLBase.EnsembleSerial())
    @test a == simulate(m, pop; seed = 7, replicates = 2, ensemble = SciMLBase.EnsembleThreads())
    @test a != simulate(m, pop; seed = 8, replicates = 2)
    # replicates differ from each other
    @test rowsof(a; id = "3", rep = 1) != rowsof(a; id = "3", rep = 2)
    # order and composition of the cohort do not change any subject's draws
    perm = shuffle(Xoshiro(1), collect(eachindex(pop)))
    b = simulate(m, pop[perm]; seed = 7, replicates = 2)
    @test b.settings == a.settings
    byid(res) = Dict((r.replicate, r.id) => r for r in res.subjects)
    ra, rb = byid(a), byid(b)
    @test all(ra[k] == rb[k] for k in keys(ra))
    for s in pop, rep in 1:2
        @test rowsof(a; id = s.id, rep) == rowsof(b; id = s.id, rep)
    end
    c = simulate(m, pop[1:10]; seed = 7, replicates = 2)                      # fewer subjects
    @test all(rowsof(a; id = s.id, rep) == rowsof(c; id = s.id, rep) for s in pop[1:10], rep in 1:2)
    # eta does not depend on whether residual error is simulated, and residual = false gives dv = ipred
    d = simulate(m, pop; seed = 7, replicates = 2, residual = false)
    @test [r.η for r in d.subjects] == [r.η for r in a.subjects]
    @test all(r -> r.dv == r.ipred, d.obs)
    @test [r.ipred for r in d.obs] == [r.ipred for r in a.obs]
    # single subject: same rng, same result
    s1 = pop[1]
    @test simulate(m, s1; rng = Xoshiro(3)) == simulate(m, s1; rng = Xoshiro(3))
    @test simulate(m, s1; rng = Xoshiro(3)).obs != simulate(m, s1; rng = Xoshiro(4)).obs
    # stream seeds depend on (seed, replicate, id, stream) only
    sd = QSPAtlas.stream_seed
    @test sd(1, 1, "a", :eta) == sd(1, 1, "a", :eta)
    @test length(unique([sd(1, 1, "a", :eta), sd(1, 1, "a", :residual), sd(1, 2, "a", :eta), sd(1, 1, "b", :eta), sd(2, 1, "a", :eta)])) == 5
    # virtual populations: covariates for subject i do not change when n grows
    cs = rng -> (WT = 60 + 20 * rand(rng),)
    v1 = virtual_population(5; regimen = dose100(), obstimes = [1.0], covariate_sampler = cs, seed = 4)
    v2 = virtual_population(8; regimen = dose100(), obstimes = [1.0], covariate_sampler = cs, seed = 4)
    @test [s.covariates for s in v1] == [s.covariates for s in v2[1:5]] && [s.id for s in v2] == string.(1:8)
end

@testset "time-varying covariate" begin
    # WT 70 -> 140 at t = 10 raises CL by 2^0.75 from then on
    m = onecmt_model()
    CL1, CL2 = 2.0, 2.0 * 2.0^0.75
    s = Subject("w"; regimen = dose100(), covariates = (WT = 70.0,), timevarying = [10.0 => (WT = 140.0,)],
        obstimes = [0.0, 5.0, 10.0, 12.0, 24.0])
    r = zero_re_prediction(m, s)
    c(tt) = tt <= 10 ? 10 * exp(-CL1 / 10 * tt) : 10 * exp(-CL1) * exp(-CL2 / 10 * (tt - 10))
    @test [o.ipred for o in r.obs] ≈ c.(s.obstimes) rtol = 1e-7
    # without the change the curve is different (the test has teeth)
    s0 = Subject("w0"; regimen = dose100(), covariates = (WT = 70.0,), obstimes = s.obstimes)
    @test zero_re_prediction(m, s0).obs[end].ipred > 1.5 * r.obs[end].ipred
    # change coinciding with a dose, and an observation at that time sees the state after both
    s2 = Subject("w2"; regimen = dose100(ii = 10.0, addl = 1), covariates = (WT = 70.0,),
        timevarying = [10.0 => (WT = 140.0,)], obstimes = [10.0, 15.0])
    r2 = zero_re_prediction(m, s2)
    top = 10 * exp(-CL1) + 10
    @test [o.ipred for o in r2.obs] ≈ [top, top * exp(-CL2 / 10 * 5)] rtol = 1e-7
    # several changes, given out of order, and one at the start folded into the baseline
    s3 = Subject("w3"; regimen = dose100(), covariates = (WT = 70.0,), obstimes = [24.0],
        timevarying = [12.0 => (WT = 70.0,), 0.0 => (WT = 140.0,), 6.0 => (WT = 35.0,)])
    @test s3.covariate_times == [0.0, 6.0, 12.0]
    CLa, CLb = 2.0 * 2.0^0.75, 2.0 * 0.5^0.75
    want = 10 * exp(-CLa / 10 * 6) * exp(-CLb / 10 * 6) * exp(-2.0 / 10 * 12)
    @test zero_re_prediction(m, s3).obs[1].ipred ≈ want rtol = 1e-7
    @test QSPAtlas.covariates_at(s3, 7.0) == (WT = 35.0,) && QSPAtlas.covariates_at(s3, -1.0) == (WT = 70.0,)
    # a change after the end of the simulation is ignored, and one with random effects still works
    s4 = Subject("w4"; regimen = dose100(), covariates = (WT = 70.0,), timevarying = [99.0 => (WT = 140.0,)], obstimes = [24.0])
    @test zero_re_prediction(m, s4).obs[1].ipred ≈ 10 * exp(-4.8) rtol = 1e-7
    r5 = simulate(m, s; η = (CL = 0.5,), residual = false)
    ηCL = exp(0.5)
    c5(tt) = tt <= 10 ? 10 * exp(-CL1 * ηCL / 10 * tt) : 10 * exp(-CL1 * ηCL) * exp(-CL2 * ηCL / 10 * (tt - 10))
    @test [o.ipred for o in r5.obs] ≈ c5.(s.obstimes) rtol = 1e-7
    # the params function must keep returning the same parameters
    @test_throws ArgumentError Subject("x"; regimen = dose100(), obstimes = [1.0], timevarying = [1.0 => (WT = 1.0,), 1.0 => (WT = 2.0,)])
end

@testset "replicates, row counts, failures" begin
    sys = onecmt
    m = PopulationModel(sys; routes = [Route(:central, sys.A)], θ = (CL = 2.0, V = 10.0), ηnames = [:CL],
        Ω = [0.09;;], params = (θ, η, c) -> [sys.CL => θ.CL * exp(η.CL), sys.V => θ.V],
        observed = (C = sys.A / sys.V, amount = sys.A), error = (C = ProportionalError(0.1), amount = AdditiveError(0.5)),
        tspan = (0.0, 24.0), solver = Tsit5(), solve_kwargs = TOL)
    pop = vcat(virtual_population(5; regimen = dose100(), obstimes = [1.0, 2.0, 2.0, 8.0]),
        [Subject("none"; regimen = dose100(), obstimes = Float64[])])
    res = simulate(m, pop; replicates = 3, seed = 1)
    @test length(res.obs) == 3 * 5 * 2 * 4                      # replicates * subjects * endpoints * times
    @test length(res.subjects) == 3 * 6
    @test sort(unique(r.replicate for r in res.obs)) == [1, 2, 3]
    @test Set(r.endpoint for r in res.obs) == Set([:C, :amount])
    @test all(r -> r.retcode === :Success, res.subjects) && all(r -> r.retcode === :Success, res.obs)
    @test count(r -> r.id == "none", res.subjects) == 3 && !any(r -> r.id == "none", res.obs)
    # duplicate observation times give duplicate rows with their own residual draws
    two = [r for r in res.obs if r.id == "1" && r.replicate == 1 && r.endpoint == :C && r.time == 2.0]
    @test length(two) == 2 && two[1].ipred == two[2].ipred && two[1].dv != two[2].dv
    # the amount and concentration observables are consistent
    amt = [r.ipred for r in res.obs if r.id == "2" && r.replicate == 2 && r.endpoint == :amount]
    con = [r.ipred for r in res.obs if r.id == "2" && r.replicate == 2 && r.endpoint == :C]
    @test amt ≈ 10 .* con
    @test res.subjects[1] == (replicate = 1, id = "1", η = res.subjects[1].η, covariates = (;), retcode = :Success)
    @test keys(res.obs[1]) == (:replicate, :id, :time, :endpoint, :ipred, :dv, :retcode)
    @test res.settings.seed == 1 && res.settings.replicates == 3 && res.settings.solver == "Tsit5"
    @test res.settings.abstol == 1e-10 && haskey(res.settings.versions, :julia)
    @test simulate(m, Subject[]).obs == NamedTuple[]
    @test occursin("SimResult", sprint(show, res))

    # a solver that gives up is reported, never dropped
    bad = PopulationModel(sys; routes = [Route(:central, sys.A)], θ = (CL = 2.0, V = 10.0), ηnames = Symbol[],
        Ω = zeros(0, 0), params = (θ, η, c) -> [sys.CL => θ.CL, sys.V => θ.V],
        observed = (C = sys.A / sys.V,), error = (C = AdditiveError(0.1),),
        tspan = (0.0, 24.0), solver = Tsit5(), solve_kwargs = (abstol = 1e-12, reltol = 1e-12, maxiters = 3))
    s = Subject("hard"; regimen = dose100(), obstimes = [0.5, 12.0, 24.0])
    r = @test_logs (:warn, r"did not finish") match_mode = :any simulate(bad, [s]; η = :sample)
    @test length(r.obs) == 3 && r.subjects[1].retcode !== :Success
    @test all(o -> o.retcode === r.subjects[1].retcode, r.obs)
    @test isnan(r.obs[end].ipred) && isnan(r.obs[end].dv)
    @test r == (@test_logs (:warn, r"did not finish") match_mode = :any simulate(bad, [s]; η = :sample))     # NaN-safe equality
end

@testset "validation" begin
    ok = (; routes = [Route(:central, onecmt.A)], θ = (CL = 2.0, V = 10.0), params = (θ, η, c) -> [onecmt.CL => θ.CL, onecmt.V => θ.V],
        observed = (C = onecmt.A / onecmt.V,), error = (C = AdditiveError(1.0),), tspan = (0.0, 24.0), solver = Tsit5())
    mk(; kw...) = PopulationModel(onecmt; merge(ok, (; ηnames = [:CL], Ω = [0.1;;]), (; kw...))...)
    @test mk() isa PopulationModel
    @test mk(ηnames = [:CL, :V], Ω = [0.1 0.0; 0.0 0.0]) isa PopulationModel        # zero variance allowed
    @test mk(ηnames = [:CL, :V], Ω = [1.0 1.0; 1.0 1.0]) isa PopulationModel        # singular PSD allowed
    @test_throws ArgumentError mk(Ω = [0.1 0.0; 0.0 0.1])                           # size mismatch
    @test_throws ArgumentError mk(ηnames = [:CL, :V], Ω = [0.1;;])
    @test_throws ArgumentError mk(ηnames = [:CL, :V], Ω = [0.1 0.05; 0.0 0.1])      # not symmetric
    @test_throws ArgumentError mk(ηnames = [:CL, :V], Ω = [1.0 2.0; 2.0 1.0])       # not PSD
    @test_throws ArgumentError mk(Ω = [-0.1;;])                                     # negative variance
    # tolerances are relative to the scale of Ω: a small indefinite matrix is still refused
    @test_throws ArgumentError mk(ηnames = [:CL, :V], Ω = [1e-12 2e-12; 2e-12 1e-12])
    @test_throws ArgumentError mk(ηnames = [:CL, :V], Ω = [1e6 2e6; 2e6 1e6])
    @test mk(ηnames = [:CL, :V], Ω = [1e-12 1e-12; 1e-12 1e-12]) isa PopulationModel   # tiny but singular PSD is fine
    @test_throws ArgumentError mk(Ω = [NaN;;])
    @test_throws ArgumentError mk(ηnames = [:CL, :CL], Ω = [0.1 0.0; 0.0 0.1])      # repeated name
    @test_throws ArgumentError mk(error = (C = AdditiveError(1.0), X = AdditiveError(1.0)))   # error key without observable
    @test_throws ArgumentError mk(error = (D = AdditiveError(1.0),))                # name mismatch
    @test_throws ArgumentError mk(error = (C = 0.1,))                               # not an ErrorModel
    @test_throws ArgumentError mk(tspan = (5.0, 1.0))
    @test_throws ArgumentError mk(routes = [Route(:central, onecmt.A), Route(:central, onecmt.A)])
    @test_throws ArgumentError Subject("s"; regimen = dose100(), obstimes = [1.0], covariates = [1.0])
    @test_throws ArgumentError Subject("s"; regimen = dose100(), obstimes = [NaN])

    m = mk()
    s = Subject(7; regimen = dose100(), obstimes = [3.0, 1.0])
    @test s.id == "7" && s.obstimes == [1.0, 3.0]
    @test_throws ArgumentError simulate(m, s; η = (V = 0.1,))                       # wrong eta names
    @test_throws ArgumentError simulate(m, s; η = [0.1, 0.2])
    @test_throws ArgumentError simulate(m, s; η = :nonsense)
    @test_throws ArgumentError simulate(m, Subject("late"; regimen = dose100(), obstimes = [30.0]); tspan = (0.0, 24.0))
    @test_throws ArgumentError simulate(m, [s, s])                                  # duplicate ids
    @test_throws ArgumentError simulate(m, [s]; replicates = 0)
    # a longer observation time extends the default tspan instead of failing
    @test simulate(m, Subject("long"; regimen = dose100(), obstimes = [30.0]); η = :zero).obs[1].retcode === :Success
    # params must return the same parameters for every subject
    flip = PopulationModel(onecmt; ok..., ηnames = [:CL], Ω = [0.1;;],
        params = (θ, η, c) -> haskey(c, :alt) ? [onecmt.V => 10.0] : [onecmt.CL => 2.0, onecmt.V => 10.0], solve_kwargs = TOL)
    @test simulate(flip, Subject("p"; regimen = dose100(), obstimes = [1.0])).obs[1].ipred > 0
    @test_throws ArgumentError simulate(flip, [Subject("p"; regimen = dose100(), obstimes = [1.0]), Subject("q"; regimen = dose100(), obstimes = [1.0], covariates = (alt = 1,))])
end

@testset "infusion, lag and bioavailability reach the model; subjects do not interfere" begin
    @parameters CLi = 2.0 Vi = 10.0 Rin = 0.0
    @variables Ai(t) = 0.0
    @mtkcompile inf = System([D(Ai) ~ -CLi / Vi * Ai + Rin], t)
    m = PopulationModel(inf; routes = [Route(:central, inf.Ai; rate = inf.Rin)], θ = (CL = 2.0, V = 10.0),
        ηnames = [:CL], Ω = [0.04;;], params = (θ, η, c) -> [inf.CLi => θ.CL * exp(η.CL), inf.Vi => θ.V],
        observed = (A = inf.Ai,), error = (A = AdditiveError(0.1),), tspan = (0.0, 24.0), solver = Tsit5(), solve_kwargs = TOL)
    times = [2.0, 5.0, 8.0, 24.0]
    k = 0.2
    A5 = 100 * (1 - exp(-k * 5))
    exact_inf(tt) = tt <= 5 ? 100 * (1 - exp(-k * tt)) : A5 * exp(-k * (tt - 5))
    s_inf = Subject("inf"; regimen = dose100(rate = 20.0), obstimes = times)             # 100 over 5 time units
    s_lag = Subject("lag"; regimen = dose100(), obstimes = times,
        control = DoseControl(bioav = Dict(:central => 0.5), lag = Dict(:central => 2.0)))
    # half the bioavailability of a fixed-rate infusion shortens it to 2.5 time units
    s_f = Subject("f"; regimen = dose100(rate = 20.0), obstimes = times, control = DoseControl(bioav = Dict(:central => 0.5)))
    exact_f(tt) = tt <= 2.5 ? 100 * (1 - exp(-k * tt)) : 100 * (1 - exp(-k * 2.5)) * exp(-k * (tt - 2.5))
    ipred(res, id) = [r.ipred for r in res.obs if r.id == id]
    pop = [s_inf, s_lag, s_f]
    res = simulate(m, pop; η = :zero, residual = false)
    @test ipred(res, "inf") ≈ exact_inf.(times) rtol = 1e-7
    @test ipred(res, "lag") ≈ [tt < 2 ? 0.0 : 50 * exp(-k * (tt - 2)) for tt in times] rtol = 1e-7 atol = 1e-12
    @test ipred(res, "f") ≈ exact_f.(times) rtol = 1e-7
    # the cohort result is the same as simulating each subject alone, and the model is reusable
    for s in pop
        @test ipred(simulate(m, s; η = :zero, residual = false), s.id) == ipred(res, s.id)
    end
    # many identical subjects solved in parallel do not share state: all give the same curve
    many = [Subject(string(i); regimen = dose100(rate = 20.0), obstimes = times) for i in 1:64]
    r = simulate(m, many; η = :zero, residual = false)
    @test all(ipred(r, string(i)) == ipred(r, "1") for i in 1:64)
    @test ipred(r, "1") ≈ exact_inf.(times) rtol = 1e-7
    @test r == simulate(m, many; η = :zero, residual = false, ensemble = SciMLBase.EnsembleSerial())
end

@testset "subject-specific initial conditions, reordered cohorts, model reuse" begin
    @parameters kic = 0.2 Vic = 10.0
    @variables Aic(t) = 1.0 Bic(t) = 0.5 * Vic            # Bic: default that depends on a parameter
    @mtkcompile icsys = System([D(Aic) ~ -kic * Aic, D(Bic) ~ 0], t)
    m = PopulationModel(icsys; routes = Route[], θ = (k = 0.2, V = 10.0), ηnames = [:k], Ω = [0.09;;],
        params = (θ, η, c) -> [icsys.kic => θ.k * exp(η.k), icsys.Vic => θ.V * c.WT / 70],
        init = (θ, η, c) -> [icsys.Aic => c.A0],            # baseline from a covariate
        observed = (A = icsys.Aic, B = icsys.Bic), error = (A = AdditiveError(0.1), B = AdditiveError(0.1)),
        tspan = (0.0, 10.0), solver = Tsit5(), solve_kwargs = TOL)
    obstimes = [0.0, 1.0, 3.0]
    none = Regimen(Dose[])
    pop = [Subject("s$i"; regimen = none, obstimes, covariates = (WT = wt, A0 = a0))
           for (i, (wt, a0)) in enumerate([(50.0, 1.0), (70.0, 5.0), (90.0, 9.0), (110.0, 2.0)])]
    function check(res)
        for sr in res.subjects
            s = pop[findfirst(p -> p.id == sr.id, pop)]
            k = 0.2 * exp(sr.η.k)
            for o in res.obs
                (o.id == sr.id && o.replicate == sr.replicate) || continue
                want = o.endpoint === :A ? s.covariates.A0 * exp(-k * o.time) : 0.5 * 10.0 * s.covariates.WT / 70
                @test o.ipred ≈ want rtol = 1e-7
            end
        end
    end
    a = simulate(m, pop; seed = 3, replicates = 2); check(a)
    b = simulate(m, reverse(pop); seed = 3, replicates = 2); check(b)             # reordered cohort
    for s in pop                                                                   # same model, one subject at a time
        r = simulate(m, s; η = :zero, residual = false)
        @test [o.ipred for o in r.obs if o.endpoint === :A] ≈ s.covariates.A0 .* exp.(-0.2 .* obstimes) rtol = 1e-7
        @test [o.ipred for o in r.obs if o.endpoint === :B] ≈ fill(5.0 * s.covariates.WT / 70, 3) rtol = 1e-7
    end
    check(simulate(m, pop[[3, 1, 4, 2]]; seed = 3))                                # and again after all that reuse
    @test a == simulate(m, pop; seed = 3, replicates = 2)
end

@testset "covariate changes leave infusions intact; route rates are not parameters" begin
    @parameters CLr = 2.0 Vr = 10.0 Rr = 0.0
    @variables Ar(t) = 0.0
    @mtkcompile infc = System([D(Ar) ~ -CLr / Vr * Ar + Rr], t)
    mk(params) = PopulationModel(infc; routes = [Route(:central, infc.Ar; rate = infc.Rr)], θ = (CL = 2.0, V = 10.0),
        ηnames = Symbol[], Ω = zeros(0, 0), params, observed = (A = infc.Ar,), error = (A = AdditiveError(0.1),),
        tspan = (0.0, 24.0), solver = Tsit5(), solve_kwargs = TOL)
    m = mk((θ, η, c) -> [infc.CLr => θ.CL * (get(c, :WT, 70.0) / 70)^0.75, infc.Vr => θ.V])
    k1, k2 = 0.2, 0.2 * 2.0^0.75
    R = 20.0
    times = [1.0, 2.0, 3.5, 5.0, 8.0]
    # change at t = 2, in the middle of a 5-unit infusion
    A2 = 100 * (1 - exp(-k1 * 2))
    A5 = R / k2 + (A2 - R / k2) * exp(-k2 * 3)
    exact(tt) = tt <= 2 ? 100 * (1 - exp(-k1 * tt)) : tt <= 5 ? R / k2 + (A2 - R / k2) * exp(-k2 * (tt - 2)) : A5 * exp(-k2 * (tt - 5))
    s = Subject("mid"; regimen = Regimen(100.0; cmt = :central, rate = R), covariates = (WT = 70.0,), obstimes = times,
        timevarying = [2.0 => (WT = 140.0,)])
    @test [o.ipred for o in zero_re_prediction(m, s).obs] ≈ exact.(times) rtol = 1e-7
    # change coinciding with the end of the infusion, and with its start
    A5b = 100 * (1 - exp(-k1 * 5))
    s2 = Subject("stop"; regimen = Regimen(100.0; cmt = :central, rate = R), covariates = (WT = 70.0,), obstimes = [4.0, 5.0, 8.0],
        timevarying = [5.0 => (WT = 140.0,)])
    @test [o.ipred for o in zero_re_prediction(m, s2).obs] ≈ [100 * (1 - exp(-k1 * 4)), A5b, A5b * exp(-k2 * 3)] rtol = 1e-7
    s3 = Subject("start"; regimen = Regimen(100.0; cmt = :central, rate = R, time = 2.0), covariates = (WT = 70.0,), obstimes = [3.0, 7.0, 8.0],
        timevarying = [2.0 => (WT = 140.0,)])
    ex3(tt) = tt <= 7 ? R / k2 * (1 - exp(-k2 * (tt - 2))) : R / k2 * (1 - exp(-k2 * 5)) * exp(-k2 * (tt - 7))
    @test [o.ipred for o in zero_re_prediction(m, s3).obs] ≈ ex3.([3.0, 7.0, 8.0]) rtol = 1e-7
    # a params function that returns the route's rate parameter is refused
    @test_throws ArgumentError mk((θ, η, c) -> [infc.CLr => θ.CL, infc.Rr => 0.0])
    # also when it needs covariates and so can only be checked at the first simulation
    late = mk((θ, η, c) -> [infc.CLr => θ.CL * c.WT / 70, infc.Rr => 0.0])
    @test_throws ArgumentError simulate(late, Subject("x"; regimen = dose100(), covariates = (WT = 70.0,), obstimes = [1.0]))
end

@testset "empty observation schedule behaves the same alone and in a cohort" begin
    m = onecmt_model()
    s = Subject("none"; regimen = dose100(), obstimes = Float64[])
    one = simulate(m, s)
    many = simulate(m, [s])
    @test isempty(one.obs) && isempty(many.obs)
    @test one.subjects[1].retcode === :Success && many.subjects[1].retcode === :Success
    # a solve that fails is reported as failed even without observation times
    bad = onecmt_model(; solve_kwargs = (abstol = 1e-12, reltol = 1e-12, maxiters = 3))
    r1 = @test_logs (:warn, r"did not finish") match_mode = :any simulate(bad, s)
    r2 = @test_logs (:warn, r"did not finish") match_mode = :any simulate(bad, [s])
    @test r1.subjects[1].retcode !== :Success && r1.subjects[1].retcode === r2.subjects[1].retcode
end

@testset "results do not depend on the number of Julia threads (separate processes)" begin
    script = joinpath(mktempdir(), "run.jl")
    write(script, """
        using QSPAtlas, ModelingToolkit, OrdinaryDiffEq
        using ModelingToolkit: t_nounits as t, D_nounits as D
        @parameters CL = 1.0 V = 10.0
        @variables A(t) = 0.0
        @mtkcompile sys = System([D(A) ~ -CL / V * A], t)
        m = PopulationModel(sys; routes = [Route(:central, sys.A)], θ = (CL = 2.0, V = 10.0), ηnames = [:CL, :V],
            Ω = [0.09 0.02; 0.02 0.04],
            params = (θ, η, c) -> [sys.CL => θ.CL * exp(η.CL), sys.V => θ.V * exp(η.V)],
            observed = (C = sys.A / sys.V,), error = (C = CombinedError(0.05, 0.1),),
            tspan = (0.0, 24.0), solver = Tsit5(), solve_kwargs = (abstol = 1e-10, reltol = 1e-10))
        pop = virtual_population(80; regimen = Regimen(100.0; cmt = :central, ii = 6.0, addl = 2), obstimes = [1.0, 4.0, 12.0, 24.0])
        res = simulate(m, pop; seed = 99, replicates = 3)
        open(ARGS[1], "w") do io
            println(io, "threads ", Threads.nthreads())
            for o in res.obs
                println(io, o.replicate, ",", o.id, ",", o.time, ",", o.endpoint, ",", o.ipred, ",", o.dv, ",", o.retcode)
            end
            for r in res.subjects
                println(io, r.replicate, ",", r.id, ",", r.η, ",", r.retcode)
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
    @test length(l1) > 1000 && l1[2:end] == l4[2:end]
end

end # module PopulationTests
