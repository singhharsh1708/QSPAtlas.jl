# Tests for Part 4: percentiles, VPC, simulation plots and the :band reference kind.
# Own module, so helper names cannot clash with other test files. The VPC model is the
# one-compartment IV bolus with a log-normal clearance, whose quantiles are known in closed form.
module SummariesTests

using QSPAtlas, Test
using ModelingToolkit, OrdinaryDiffEq, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D
using Random: Xoshiro, randn
using Statistics: median, quantile, std, mean

# ---- helpers on hand-built results -----------------------------------------------------------

row(rep, id, time, v; endpoint = :C) = (replicate = rep, id = string(id), time = Float64(time), endpoint, ipred = Float64(v), dv = Float64(v), retcode = :Success)
simresult(obs) = SimResult(obs, [(replicate = r, id = "s", η = (;), covariates = (;), retcode = :Success) for r in unique(o.replicate for o in obs)], (;))
mean_(x) = sum(x) / length(x)
at(rows; kw...) = only([r for r in rows if all(getproperty(r, k) == v for (k, v) in kw)])

@testset "percentiles: hand-computed values" begin
    # three times, five subjects, one replicate
    cols = Dict(1 => [1.0, 2, 3, 4, 10], 2 => [5.0, 6, 7, 8, 9], 3 => [11.0, 12, 13, 14, 15])
    r = simresult([row(1, s, tm, v) for (tm, vs) in cols for (s, v) in enumerate(vs)])
    p = percentiles(r, :C; q = (0.1, 0.5, 0.9))
    @test length(p) == 9
    # Hyndman-Fan type 7: h = (n - 1) q + 1 on the sorted values
    @test at(p; t_mid = 1.0, q = 0.1).value ≈ 1.4              # 1 + 0.4 * (2 - 1)
    @test at(p; t_mid = 1.0, q = 0.5).value ≈ 3.0
    @test at(p; t_mid = 1.0, q = 0.9).value ≈ 7.6              # 4 + 0.6 * (10 - 4)
    @test at(p; t_mid = 2.0, q = 0.9).value ≈ 8.6
    @test all(x -> x.n == 5 && x.n_failed == 0 && x.bin_lo == x.bin_hi == x.t_mid, p)
    @test [x.bin for x in p] == [1, 1, 1, 2, 2, 2, 3, 3, 3]
    # value = :ipred reads the other column
    r2 = simresult([merge(o, (; dv = 100 * o.dv)) for o in r.obs])
    @test at(percentiles(r2, :C; q = (0.5,)); t_mid = 1.0).value == 300.0
    @test at(percentiles(r2, :C; q = (0.5,), value = :ipred); t_mid = 1.0).value == 3.0

    # edges: [lo, hi) with the last bin closed, so an observation on an inner edge goes up
    pe = percentiles(r, :C; q = (0.5,), bins = [1.0, 2.0, 3.0])
    @test [(x.bin_lo, x.bin_hi, x.n) for x in pe] == [(1.0, 2.0, 5), (2.0, 3.0, 10)]
    @test pe[2].value ≈ 10.0                                    # t = 2 and t = 3 pooled
    @test pe[1].t_mid == 1.5
    @test only(percentiles(r, :C; q = (0.5,), bins = [0.0, 3.0])).value ≈ 8.0     # all 15 values
    @test only(percentiles(r, :C; q = (0.5,), bins = [0.0, 1.0])).n == 5          # last edge included
    @test only(percentiles(r, :C; q = (0.5,), bins = [1.5, 2.5])).n == 5          # only t = 2
    @test only(percentiles(r, :C; q = (0.5,), bins = [3.5, 4.5])).n == 0          # empty bin
    @test isnan(only(percentiles(r, :C; q = (0.5,), bins = [3.5, 4.5])).value)

    # per replicate
    r3 = simresult([row(rep, s, 1.0, v) for (rep, off) in ((1, 0.0), (2, 10.0)) for (s, v) in enumerate([1.0, 2, 3, 4, 5] .+ off)])
    pr = percentiles(r3, :C; q = (0.5,), by_replicate = true)
    @test [(x.replicate, x.value) for x in pr] == [(1, 3.0), (2, 13.0)]
    @test percentiles(r3, :C; q = (0.5,))[1].value ≈ 8.0 && percentiles(r3, :C; q = (0.5,))[1].n == 10

    # a second endpoint is not mixed in
    both = simresult(vcat(r.obs, [row(1, s, 1.0, 1000.0 + s; endpoint = :E) for s in 1:5]))
    @test at(percentiles(both, :C; q = (0.5,)); t_mid = 1.0).value == 3.0
    @test at(percentiles(both, :E; q = (0.5,))).value == 1003.0

    # bad input
    @test_throws ArgumentError percentiles(r, :nope)
    @test_throws ArgumentError percentiles(r, :C; value = :foo)
    @test_throws ArgumentError percentiles(r, :C; q = (0.5, 1.2))
    @test_throws ArgumentError percentiles(r, :C; q = (0.5, 0.5))
    @test_throws ArgumentError percentiles(r, :C; bins = [1.0, 1.0, 2.0])
    @test_throws ArgumentError percentiles(r, :C; bins = [1.0])
    @test_throws ArgumentError percentiles(r, :C; bins = :something)
    @test_throws ArgumentError percentiles(r, :C; blq = (lloq = 1.0, handling = :magic))
end

@testset "percentiles: failed and NaN values are excluded and counted" begin
    vs = [1.0, 2.0, NaN, 4.0, 10.0]
    r = simresult([row(1, s, 1.0, v) for (s, v) in enumerate(vs)])
    p = @test_logs (:warn, r"non-finite") percentiles(r, :C; q = (0.5,))
    @test p[1].value ≈ 3.0 && p[1].n == 4 && p[1].n_failed == 1       # median of [1, 2, 4, 10]
    # too few usable values: NaN quantile, n still reported
    p = @test_logs (:warn, r"non-finite") percentiles(r, :C; q = (0.5,), min_n = 5)
    @test isnan(p[1].value) && p[1].n == 4 && p[1].n_failed == 1
    # an all-failed bin
    r = simresult([row(1, s, 1.0, NaN) for s in 1:3])
    p = @test_logs (:warn, r"non-finite") percentiles(r, :C; q = (0.1, 0.9))
    @test all(x -> isnan(x.value) && x.n == 0 && x.n_failed == 3, p)
    # infinities count as failed too
    r = simresult([row(1, s, 1.0, v) for (s, v) in enumerate([1.0, Inf, 3.0])])
    p = @test_logs (:warn, r"non-finite") percentiles(r, :C; q = (0.5,))
    @test p[1].value ≈ 2.0 && p[1].n_failed == 1
end

@testset "percentiles: BLQ handling" begin
    r = simresult([row(1, s, 1.0, v) for (s, v) in enumerate([0.5, 0.8, 2.0, 3.0, 4.0])])
    f(h; q = 0.25) = only(percentiles(r, :C; q = (q,), blq = h === nothing ? nothing : (lloq = 1.0, handling = h)))
    @test f(nothing).value ≈ 0.8 && f(nothing).n_blq == 0
    k = f(:keep)
    @test k.value ≈ 0.8 && k.n == 5 && k.n_blq == 2
    c = f(:censor_to_lloq)
    @test c.value ≈ 1.0 && c.n == 5 && c.n_blq == 2                 # [1, 1, 2, 3, 4]
    d = f(:drop)
    @test d.value ≈ 2.5 && d.n == 3 && d.n_blq == 2                 # [2, 3, 4]
    # default handling is censoring; a value exactly at lloq is not BLQ
    @test only(percentiles(r, :C; q = (0.25,), blq = (lloq = 1.0,))).value ≈ 1.0
    @test only(percentiles(r, :C; q = (0.5,), blq = (lloq = 2.0, handling = :drop))).n_blq == 2
    @test only(percentiles(r, :C; q = (0.5,), blq = (lloq = 2.0, handling = :drop))).n == 3
end

# ---- the MC standard error -------------------------------------------------------------------

@testset "mc_se_median (bootstrap over replicates)" begin
    rng = Xoshiro(11)
    se(x; seed = 1) = QSPAtlas.mc_se_median(x, Xoshiro(seed))
    @test isnan(se(ones(9)))
    @test se(ones(50)) == 0.0                                  # no spread at all: really no error
    x = randn(rng, 300)
    @test se(x; seed = 5) == se(x; seed = 5)                   # reproducible from the generator
    @test se(x; seed = 5) != se(x; seed = 6)
    # normal data: close to the textbook 1.2533 / sqrt(K), on average
    for K in (100, 400)
        @test mean_([se(randn(rng, K); seed = i) for i in 1:100]) ≈ 1.2533 / sqrt(K) rtol = 0.12
    end
    # discrete values (BLQ censoring ties): 61 zeros and 39 ones. A rank interval has width 0 here
    # and would report no error; the true sd of the median of such draws is about 0.14.
    @test se([zeros(61); ones(39)]) > 0.05
end

# independent batches: the sd of the across-replicate median over independent vpc runs (different
# seeds, fresh replicates) against the mean reported mc_se. Built from random draws directly, so it is cheap.
@testset "mc_se against independent batches" begin
    function batch(b; K = 100, discrete = false)
        rng = Xoshiro(1000 + b)
        obs = [row(k, 1, 1.0, discrete ? Float64(rand(rng) > 0.5) : exp(0.3 * randn(rng))) for k in 1:K]
        return only(vpc(simresult(obs); endpoint = :C, q = (0.5,), seed = b).rows)
    end
    for discrete in (false, true)
        xs = [batch(b; discrete) for b in 1:300]
        sd = std([x.median for x in xs])
        # a bootstrap that cannot estimate the SE (heavily tied) now gives NaN instead of 0 and is left
        # out of the mean; it must stay rare here (a balanced 0/1 sample, 100 replicates)
        ses = [x.mc_se for x in xs]
        @test count(isnan, ses) <= 0.03 * length(ses)
        reported = mean_(filter(isfinite, ses))
        @test sd / 1.5 < reported < sd * 1.5
    end
end

# ---- vpc on hand-built results ---------------------------------------------------------------

@testset "vpc: summary of replicate quantiles" begin
    # replicate k has three subjects k - 1, k, k + 1, so its median is k
    K = 20
    r = simresult([row(k, s, 1.0, k + off) for k in 1:K for (s, off) in enumerate((-1, 0, 1))])
    v = vpc(r; endpoint = :C, q = (0.5,))
    x = only(v.rows)
    @test x.median ≈ 10.5 && x.n == 3 && x.n_rep == K && x.n_failed == 0
    @test x.lo ≈ 1.475 && x.hi ≈ 19.525                             # type 7 at 0.025 and 0.975 of 1:20
    @test 0.5 * 1.2533 * std(1:K) / sqrt(K) < x.mc_se < 2 * 1.2533 * std(1:K) / sqrt(K)   # bootstrap, normal-theory scale
    @test v.replicates == K && v.observed === nothing && v.pi_level == 0.95
    v90 = vpc(r; endpoint = :C, q = (0.5,), pi_level = 0.5)
    @test only(v90.rows).lo ≈ quantile(1:K, 0.25) && only(v90.rows).hi ≈ quantile(1:K, 0.75)
    @test length(v) == 1 && v[1] === x && collect(v) == v.rows
    @test occursin("VPCResult", sprint(show, v))
    # two bins by edges and three quantiles keep the (bin, q) ordering
    r2 = simresult(vcat(r.obs, [row(k, s, 2.0, 100 + k + off) for k in 1:K for (s, off) in enumerate((-1, 0, 1))]))
    v2 = vpc(r2; endpoint = :C, bins = [0.5, 1.5, 2.5])
    @test [(x.bin, x.q) for x in v2.rows] == [(b, q) for b in 1:2 for q in (0.1, 0.5, 0.9)]
    @test at(v2.rows; bin = 2, q = 0.5).median ≈ 110.5

    # replicate 3 is entirely failed: excluded, counted, reported
    bad = simresult([row(k, s, 1.0, k == 3 ? NaN : k + off) for k in 1:K for (s, off) in enumerate((-1, 0, 1))])
    vb = @test_logs (:warn, r"non-finite") vpc(bad; endpoint = :C, q = (0.5,))
    xb = only(vb.rows)
    @test xb.n_rep == K - 1 && xb.n_failed == 3 && vb.n_failed == 3
    @test xb.median ≈ median(filter(!=(3), 1:K))
    # one NaN subject in one replicate: still defined, one value less
    one = simresult([row(k, s, 1.0, (k == 3 && s == 2) ? NaN : k + off) for k in 1:K for (s, off) in enumerate((-1, 0, 1))])
    vo = @test_logs (:warn, r"non-finite") vpc(one; endpoint = :C, q = (0.5,))
    @test only(vo.rows).n_rep == K && only(vo.rows).n_failed == 1
    # min_n makes the thin replicate undefined
    vm = @test_logs (:warn, r"non-finite") vpc(one; endpoint = :C, q = (0.5,), min_n = 3)
    @test only(vm.rows).n_rep == K - 1
    # too few replicates for an MC standard error
    @test isnan(only(vpc(simresult([row(k, s, 1.0, k + s) for k in 1:5 for s in 1:3]); endpoint = :C, q = (0.5,)).rows).mc_se)

    # BLQ
    rb = simresult([row(k, s, 1.0, v) for k in 1:12 for (s, v) in enumerate([0.2, 0.5, 5.0, 6.0, 7.0])])
    for (h, med, n) in ((:censor_to_lloq, 5.0, 5), (:keep, 5.0, 5), (:drop, 6.0, 3))
        x = only(vpc(rb; endpoint = :C, q = (0.5,), blq = (lloq = 1.0, handling = h)).rows)
        @test x.median ≈ med && x.n == n && x.n_blq == 2 * 12
    end
    x = only(vpc(rb; endpoint = :C, q = (0.2,), blq = (lloq = 1.0, handling = :censor_to_lloq)).rows)
    @test x.median ≈ 1.0                                           # [1, 1, 5, 6, 7] at 0.2 (h = 1.8)
    x = only(vpc(rb; endpoint = :C, q = (0.2,), blq = (lloq = 1.0, handling = :keep)).rows)
    @test x.median ≈ 0.2 + 0.8 * 0.3

    # observed data: pooled quantiles per bin, diagnostic only
    obs = [(id = "a$i", time = 1.0, dv = Float64(i)) for i in 1:5]
    push!(obs, (id = "z", time = 7.0, dv = 1.0))                     # not a simulated time
    vo = vpc(r; endpoint = :C, q = (0.5,), observed = obs)
    @test only(vo.observed).value == 3.0 && only(vo.observed).n == 5 && vo.n_unbinned == 1
    @test vpc(r; endpoint = :C, q = (0.5,), observed = obs, bins = [0.0, 10.0]).n_unbinned == 0
    @test only(vpc(r; endpoint = :C, q = (0.5,), observed = [(id = "a", time = 1.0, dv = NaN), (id = "b", time = 1.0, dv = 4.0)]).observed).n_failed == 1

    # bad options are caught before any simulation
    @test_throws ArgumentError vpc(r; endpoint = :C, pi_level = 1.0)
    @test_throws ArgumentError vpc(r; endpoint = :C, value = :x)
    @test_throws ArgumentError vpc(r; endpoint = :nope)
end

# ---- vpc on a real model with closed-form quantiles --------------------------------------------

@parameters CL = 1.0 V = 10.0
@variables A(t) = 0.0
@mtkcompile onecmt = System([D(A) ~ -CL / V * A], t)

const θCL, V0, ω, DOSE = 2.0, 10.0, 0.3, 100.0

function onecmt_model(; error = ProportionalError(0.1))
    sys = onecmt
    return PopulationModel(sys;
        routes = [Route(:central, sys.A)],
        θ = (CL = θCL, V = V0), ηnames = [:CL], Ω = [ω^2;;],
        params = (θ, η, c) -> [sys.CL => θ.CL * exp(η.CL), sys.V => θ.V],
        observed = (C = sys.A / sys.V,), error = (C = error,),
        tspan = (0.0, 12.0), solver = Tsit5(), solve_kwargs = (abstol = 1e-10, reltol = 1e-10))
end

cohort(n; times = [1.0, 4.0, 12.0]) = [Subject("s$i"; regimen = Regimen(DOSE; cmt = :central), obstimes = times) for i in 1:n]

# q-quantile of ipred = (D/V) exp(-CL e^η t / V): decreasing in η, so it is the (1 - q) quantile of η
Φ(z) = 0.5 * (1 + erf_(z / sqrt(2)))
function erf_(x)                                   # Abramowitz-Stegun 7.1.26, |error| < 1.5e-7
    a = (0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429)
    u = 1 / (1 + 0.3275911 * abs(x))
    y = 1 - (((((a[5] * u + a[4]) * u) + a[3]) * u + a[2]) * u + a[1]) * u * exp(-x^2)
    return x >= 0 ? y : -y
end
function probit(p)
    lo, hi = -8.0, 8.0
    for _ in 1:80
        mid = (lo + hi) / 2
        Φ(mid) < p ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
ipred_at(tm, η) = DOSE / V0 * exp(-θCL * exp(η) * tm / V0)
exact_ipred_quantile(tm, q) = ipred_at(tm, ω * probit(1 - q))

@testset "vpc against the closed form of a log-normal clearance" begin
    m = onecmt_model()
    N, K = 50, 100
    subjects = cohort(N)
    res = simulate(m, subjects; replicates = K, seed = 7, residual = false)
    v = vpc(res; endpoint = :C, value = :ipred)
    @test length(v) == 9 && v.replicates == K
    @test probit(0.9) ≈ 1.2815515655 atol = 1e-6
    # (a) population closed form. The type 7 sample quantile of N subjects is not exactly the
    # population quantile; its median over replicates is within 1/N IN PROBABILITY of q (stated
    # bound, N = 50 so 0.02), so the simulated median must lie between the closed form at q - 1/N
    # and q + 1/N (evaluated exactly, no derivative approximation), widened by 4 Monte Carlo SE.
    for x in v.rows
        e1, e2 = exact_ipred_quantile(x.t_mid, clamp(x.q - 1 / N, 0, 1)), exact_ipred_quantile(x.t_mid, clamp(x.q + 1 / N, 0, 1))
        @test min(e1, e2) - 4 * x.mc_se <= x.median <= max(e1, e2) + 4 * x.mc_se
        @test x.lo < x.median < x.hi && x.n == N && x.n_rep == K
    end
    # (b) sharp check with the finite-N bias included: the same estimator (type 7 quantile of N
    # subjects, median over replicates) applied to directly drawn η, no ODE and no vpc code.
    R = 4000
    rng = Xoshiro(99)
    ηs = [ω .* randn(rng, N) for _ in 1:R]
    for x in v.rows
        per = [quantile([ipred_at(x.t_mid, e) for e in η], x.q) for η in ηs]
        ref, ref_se = median(per), 1.2533 * std(per) / sqrt(R)
        @test abs(x.median - ref) <= 4 * sqrt(x.mc_se^2 + ref_se^2)
    end
    # the median line is about the median of the population at every time
    for x in filter(x -> x.q == 0.5, v.rows)
        @test x.median ≈ DOSE / V0 * exp(-θCL * x.t_mid / V0) rtol = 0.02
    end

    # mc_se matches the measured spread of the replicate quantiles and shrinks like 1/sqrt(K);
    # a small cohort keeps this cheap (the SE of the median of K values does not need many subjects)
    K2 = 400
    res2 = simulate(m, cohort(10); replicates = K2, seed = 8, residual = false)
    v2 = vpc(res2; endpoint = :C, value = :ipred, q = (0.5,))
    per_rep = percentiles(res2, :C; q = (0.5,), value = :ipred, by_replicate = true)
    for tm in (1.0, 4.0, 12.0)
        vals = [x.value for x in per_rep if x.t_mid == tm]
        @test at(v2.rows; t_mid = tm).mc_se ≈ 1.2533 * std(vals) / sqrt(K2) rtol = 0.3
    end
    sub(k) = SimResult(filter(o -> o.replicate <= k, res2.obs), filter(s -> s.replicate <= k, res2.subjects), res2.settings)
    se(k) = mean_([x.mc_se / x.median for x in vpc(sub(k); endpoint = :C, value = :ipred).rows])   # relative, times differ
    @test 0.3 < se(100) / se(25) < 0.75                            # expect 0.5
    @test 0.3 < se(400) / se(100) < 0.75
    # the interval does not shrink with more replicates
    w(k) = mean_([(x.hi - x.lo) / x.median for x in vpc(sub(k); endpoint = :C, value = :ipred).rows])
    @test w(400) ≈ w(100) rtol = 0.2

    # with the residual error the spread grows, the cohort size is unchanged
    vd = vpc(m, subjects; endpoint = :C, replicates = 30, seed = 7)
    @test vd.value === :dv && all(x -> x.n == N, vd.rows)
    @test (at(vd.rows; t_mid = 4.0, q = 0.9).hi - at(vd.rows; t_mid = 4.0, q = 0.1).lo) > 0
end

@testset "vpc is reproducible from the seed" begin
    m = onecmt_model()
    s = cohort(6)
    a = vpc(m, s; endpoint = :C, replicates = 12, seed = 3)
    b = vpc(m, s; endpoint = :C, replicates = 12, seed = 3)
    c = vpc(m, s; endpoint = :C, replicates = 12, seed = 4)
    d = vpc(m, s; endpoint = :C, replicates = 12, seed = 3, ensemble = SciMLBase.EnsembleSerial())
    @test isequal(a.rows, b.rows) && a.seed == 3
    @test isequal(a.rows, d.rows)                                   # independent of the threading
    @test !isequal(a.rows, c.rows)
    # the same cohort in another order gives the same quantiles
    e = vpc(m, reverse(s); endpoint = :C, replicates = 12, seed = 3)
    @test all(x -> isapprox(x[1].median, x[2].median; rtol = 1e-12), zip(a.rows, e.rows))
    # observed data through the model method
    obs = [(id = "o$i", time = 4.0, dv = 3.0 + 0.1 * i) for i in 1:8]
    vo = vpc(m, s; endpoint = :C, replicates = 12, seed = 3, observed = obs, q = (0.5,))
    @test at(vo.observed; t_mid = 4.0).value ≈ median([o.dv for o in obs]) && vo.n_unbinned == 0
    # option errors come before the simulation
    @test_throws ArgumentError vpc(m, s; endpoint = :C, replicates = 0)
    @test_throws ArgumentError vpc(m, s; endpoint = :C, replicates = 5, q = (2.0,))
end

# ---- the :band reference kind ----------------------------------------------------------------

ax(scale = :linear, lo = 0.0, hi = 10.0) = QSPAtlas.Axis("a", "u", scale, lo, hi)
band_ref(xs, lo, hi; xaxis = ax(), yaxis = ax()) =
    ReferenceCurve("b", "band", Float64.(xs), (Float64.(lo) .+ Float64.(hi)) ./ 2, xaxis, yaxis, :band, Dict{String, String}(), Float64.(lo), Float64.(hi))

@testset "band comparison" begin
    ref = band_ref([2.0, 4.0, 6.0], [1.0, 3.0, 5.0], [3.0, 5.0, 7.0])           # y in [x - 1, x + 1]
    xs = collect(0.0:1.0:10.0)
    c = compare_curve(ref, xs, xs)                                              # y = x is inside
    @test c.d == [0.0, 0.0, 0.0] && passes(c) && c.rms == 0
    # upper edge missed by 1 / 10 of the axis (y = x + 2 is 1 above x + 1), then below the lower edge
    @test compare_curve(ref, xs, xs .+ 2).d ≈ [0.1, 0.1, 0.1]
    @test compare_curve(ref, xs, xs .- 2).d ≈ [0.1, 0.1, 0.1]
    @test !passes(compare_curve(ref, xs, xs .+ 2))
    # on the edge is inside; just outside is a small error
    @test compare_curve(ref, xs, xs .+ 1).d == [0.0, 0.0, 0.0]
    @test compare_curve(ref, xs, xs .+ 1.05).d ≈ [0.005, 0.005, 0.005]
    # mixed: only the middle reference point is outside
    ys = copy(xs); ys[5] = 4.0 + 1.5
    @test compare_curve(ref, xs, ys).d ≈ [0.0, 0.05, 0.0]
    # interpolation between simulated points (reference x = 4 lies between 3 and 5)
    @test compare_curve(ref, [0.0, 3.0, 5.0, 10.0], [0.0, 3.0, 5.0, 10.0]).d == [0.0, 0.0, 0.0]
    @test compare_curve(ref, [0.0, 3.0, 5.0, 10.0], [0.0, 3.0, 9.0, 10.0]).d[2] ≈ 0.1   # 6 at x = 4, upper edge 5
    # Interpolate first, cap the error afterwards. Sim (0, 0) -> (10, 100) on axes [0, 10]: at x = 2.5 the
    # interpolant is 25, far above the band [4.9, 5.1]; at x = 0.5 it is 5, inside. Clipping the
    # endpoints before interpolating gets both wrong.
    cb = band_ref([2.5, 0.5], [4.9, 4.9], [5.1, 5.1])
    cc = compare_curve(cb, [0.0, 10.0], [0.0, 100.0])
    @test cc.d == [1.0, 0.0] && !passes(cc)
    @test compare_curve(band_ref([2.5], [4.9], [5.1]), [0.0, 10.0], [0.0, 100.0]).d == [1.0]
    @test compare_curve(band_ref([0.5], [4.9], [5.1]), [0.0, 10.0], [0.0, 100.0]).d == [0.0]
    @test compare_curve(band_ref([0.5], [5.2], [5.4]), [0.0, 10.0], [0.0, 100.0]).d ≈ [0.02]   # 5.0 vs lower edge 5.2
    # error is capped at 1
    far = band_ref([5.0], [0.0], [1.0]; yaxis = ax(:linear, 0.0, 10.0))
    @test compare_curve(far, [0.0, 10.0], [100.0, 100.0]).d == [1.0]

    # x outside the simulated range is the maximum error
    r2 = band_ref([2.0, 12.0], [0.0, 0.0], [10.0, 10.0]; xaxis = ax(:linear, 0.0, 20.0))
    c2 = compare_curve(r2, [0.0, 5.0, 10.0], [1.0, 1.0, 1.0])
    @test c2.d == [0.0, 1.0] && !passes(c2)
    c2 = compare_curve(r2, [3.0, 5.0, 14.0], [1.0, 1.0, 1.0])
    @test c2.d == [1.0, 0.0]
    @test compare_curve(r2, [2.0, 5.0, 12.0], [1.0, 1.0, 1.0]).d == [0.0, 0.0]        # ends count

    # a gap is not bridged: a reference x inside the gap has no simulated value
    rg = band_ref([1.0, 5.0, 9.0], [0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
    cg = compare_curve(rg, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0], [1.0, 1.0, NaN, NaN, 1.0, 1.0])
    @test cg.d == [0.0, 1.0, 0.0]
    cg = compare_curve(rg, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0], [1.0, 1.0, 1.0, Inf, 1.0, 1.0])
    @test cg.d == [0.0, 1.0, 0.0]                              # 5 lies between 4 and 6, and 6 is Inf
    # one finite point is enough when it matches a reference x exactly
    @test compare_curve(rg, [0.0, 1.0, 2.0], [NaN, 1.0, NaN]).d == [0.0, 1.0, 1.0]
    @test_throws ArgumentError compare_curve(rg, [0.0, 1.0, 2.0], [NaN, NaN, NaN])
    @test_throws ArgumentError compare_curve(rg, [2.0, 1.0, 3.0], [1.0, 1.0, 1.0])  # unsorted x
    # an exact finite match needs no neighbour, with or without a valid segment elsewhere
    r5 = band_ref([5.0], [4.0], [6.0])
    @test compare_curve(r5, [0.0, 5.0, 10.0], [NaN, 5.0, NaN]).d == [0.0]
    @test compare_curve(r5, [0.0, 5.0, 10.0], [NaN, 9.0, NaN]).d ≈ [0.3]
    @test compare_curve(band_ref([5.0, 8.0], [4.0, 0.0], [6.0, 9.0]), [0.0, 5.0, 7.0, 7.5, 10.0], [NaN, 5.0, NaN, 1.0, 1.0]).d == [0.0, 0.0]
    # bridging 4 -> 10 would put the curve at 5 inside the band [4, 6] at x = 5; the gap rule does not
    rb = band_ref([5.0], [4.0], [6.0])
    @test compare_curve(rb, [0.0, 4.0, 6.0, 10.0], [0.0, 4.0, 6.0, 10.0]).d == [0.0]
    @test compare_curve(rb, [0.0, 4.0, 6.0, 10.0], [0.0, 4.0, NaN, 10.0]).d == [1.0]

    # log axes: errors in decades; interpolation is linear in the log coordinate
    ly = ax(:log, 0.01, 100.0)                                                  # 4 decades
    rl = band_ref([5.0], [1.0], [10.0]; yaxis = ly)
    @test compare_curve(rl, [0.0, 10.0], [3.0, 3.0]).d == [0.0]
    @test compare_curve(rl, [0.0, 10.0], [100.0, 100.0]).d ≈ [0.25]             # one decade over 10, in 4 decades
    @test compare_curve(rl, [0.0, 10.0], [0.01, 0.01]).d ≈ [0.5]                # two decades under 1
    # log-linear interpolation: between 0.1 at x = 0 and 10 at x = 10 the geometric midpoint is 1
    rm = band_ref([5.0], [0.9], [1.1]; yaxis = ly)
    @test compare_curve(rm, [0.0, 10.0], [0.1, 10.0]).d ≈ [0.0] atol = 1e-12
    # log x axis: the midpoint of 1 and 100 is 10
    rx = band_ref([10.0], [4.9], [5.1]; xaxis = ax(:log, 0.1, 1000.0))
    @test compare_curve(rx, [1.0, 100.0], [0.0, 10.0]).d ≈ [0.0] atol = 1e-12

    # non-positive simulated values on a log axis are a large error, not an exception
    @test compare_curve(rl, [0.0, 10.0], [0.0, 0.0]).d[1] > 0.2

    # the tolerance still decides: rms 0.005 passes, 0.1 fails
    rc = band_ref([2.0, 4.0, 6.0], [1.0, 3.0, 5.0], [3.0, 5.0, 7.0])
    @test passes(compare_curve(rc, xs, xs .+ 1.05)) && !passes(compare_curve(rc, xs, xs .+ 2))
    @test !passes(compare_curve(rc, xs, xs .+ 1.05; tolerance = Tolerance(rms = 0.001, p95 = 0.001)))
end

@testset "vertical metric for percentile lines" begin
    vref(xs, ys; meta = Dict("metric" => "vertical"), kw...) = ReferenceCurve("v", "line", Float64.(xs), Float64.(ys), ax(), ax(), :curve, meta)
    # a steep, time-shifted line: the nearest-polyline default is fooled, the vertical metric is not.
    # sim: step at x = 5 (from (5, 0) to (5.2, 10)); the reference point (5.1, 0) lies 0.1 / 10 from the
    # polyline but the simulated value at x = 5.1 is 5.
    xs, ys = [0.0, 5.0, 5.2, 10.0], [0.0, 0.0, 10.0, 10.0]
    near = vref([5.1], [0.0]; meta = Dict{String, String}())
    @test passes(compare_curve(near, xs, ys))                                   # default unchanged
    vert = vref([5.1], [0.0])
    cv = compare_curve(vert, xs, ys)
    @test cv.d ≈ [0.5] && !passes(cv)
    # same rules as bands: exact match first, interpolation, gaps, range, cap, log axes, singleton
    @test compare_curve(vref([2.0, 4.0], [2.0, 4.5]), [0.0, 2.0, 4.0, 10.0], [0.0, 2.0, 4.0, 10.0]).d ≈ [0.0, 0.05]
    @test compare_curve(vref([5.0], [5.0]), [0.0, 5.0, 10.0], [NaN, 5.0, NaN]).d == [0.0]
    @test compare_curve(vref([5.0], [5.0]), [0.0, 4.0, 6.0, 10.0], [0.0, 4.0, NaN, 10.0]).d == [1.0]   # gap
    @test compare_curve(vref([12.0], [5.0]), [0.0, 10.0], [0.0, 10.0]).d == [1.0]                    # out of range
    @test compare_curve(vref([5.0], [0.0]), [0.0, 10.0], [100.0, 100.0]).d == [1.0]                  # capped
    ly = ax(:log, 0.01, 100.0)
    rl = ReferenceCurve("v", "line", [5.0], [1.0], ax(), ly, :curve, Dict("metric" => "vertical"))
    @test compare_curve(rl, [0.0, 10.0], [0.1, 10.0]).d ≈ [0.0] atol = 1e-12                         # log midpoint is 1
    @test compare_curve(rl, [0.0, 10.0], [10.0, 10.0]).d ≈ [0.25]                                    # one decade
    # file format: header key, round trip, validation
    p = write_reference(tempname() * ".csv", vert)
    @test occursin("# metric: vertical", read(p, String)) && load_reference(p).meta["metric"] == "vertical"
    @test passes(compare_curve(load_reference(p), xs, ys)) == false
    bad = tempname() * ".csv"
    write(bad, "# id: b\n# figure: f\n# x: t [h] linear 0 1\n# y: v [u] linear 0 1\n# kind: curve\n# metric: sideways\nx,y\n0,1\n")
    @test_throws ArgumentError load_reference(bad)
    write(bad, "# id: b\n# figure: f\n# x: t [h] linear 0 1\n# y: v [u] linear 0 1\n# kind: points\n# metric: vertical\nx,y\n0,1\n")
    @test_throws ArgumentError load_reference(bad)
end

@testset "band references: file format and constructors" begin
    # the old eight-argument constructor still works and leaves ylo/yhi empty
    r = ReferenceCurve("c", "f", [0.0, 1.0], [1.0, 2.0], ax(), ax(), :curve, Dict{String, String}())
    @test isempty(r.ylo) && isempty(r.yhi)
    @test_throws ArgumentError ReferenceCurve("c", "f", [0.0], [1.0], ax(), ax(), :curve, Dict{String, String}(), [0.0], [1.0])
    @test_throws ArgumentError band_ref([1.0], [3.0], [2.0])                       # lo > hi
    @test_throws ArgumentError ReferenceCurve("c", "f", [0.0, 1.0], [1.0, 2.0], ax(), ax(), :band, Dict{String, String}(), [0.0], [1.0])

    b = band_ref([0.0, 1.0, 2.5], [1.0, 2.0, 3.0], [2.0, 4.0, 3.5]; yaxis = ax(:log, 0.1, 100.0))
    path = write_reference(tempname() * ".csv", b)
    txt = read(path, String)
    @test occursin("# kind: band", txt) && occursin("x,lo,hi", txt) && !occursin("x,y", txt)
    b2 = load_reference(path)
    @test b2.kind === :band && b2.x == b.x && b2.ylo == b.ylo && b2.yhi == b.yhi
    @test b2.y == (b.ylo .+ b.yhi) ./ 2 && b2.yaxis.scale === :log && b2.id == "b"
    # a hand-written file
    p = tempname() * ".csv"
    write(p, "# id: fig2_pk_band\n# figure: Fig 2\n# x: time [h] linear 0 24\n# y: conc [mg/L] log 0.1 100\n# kind: band\n# source: test\nx,lo,hi\n0,1,3\n12,0.5,2\n")
    h = load_reference(p)
    @test h.kind === :band && h.ylo == [1.0, 0.5] && h.yhi == [3.0, 2.0] && h.y == [2.0, 1.25] && h.meta["source"] == "test"
    # round trip of ordinary kinds is unchanged
    c = load_reference(write_reference(tempname() * ".csv", r))
    @test c.kind === :curve && c.x == r.x && c.y == r.y && isempty(c.ylo)
    # wrong column count for the kind
    bad = tempname() * ".csv"
    write(bad, "# id: b\n# figure: f\n# x: t [h] linear 0 1\n# y: v [u] linear 0 1\n# kind: band\nx,y\n0,1\n")
    @test_throws ArgumentError load_reference(bad)
    write(bad, "# id: b\n# figure: f\n# x: t [h] linear 0 1\n# y: v [u] linear 0 1\n# kind: curve\nx,lo,hi\n0,1,2\n")
    @test_throws ArgumentError load_reference(bad)
    write(bad, "# id: b\n# figure: f\n# x: t [h] linear 0 1\n# y: v [u] linear 0 1\n# kind: band\nx,lo,hi\n0,2,1\n")
    @test_throws ArgumentError load_reference(bad)

    # a band flows through run_curves, the score file and the reproduction plot
    xs = collect(0.0:0.5:24.0)
    run = run_curves(_ -> (xs, 2.0 .* ones(length(xs))), [h], "as published")
    @test passes(run.curves[1].comparison)
    rep = Reproduction(id = "bandtest", citation = "c", doi = "d", published = run)
    s = read_score(write_score(tempname() * ".toml", rep))
    @test s["published"]["curves"][1]["passed"] && s["published"]["curves"][1]["n_points"] == 2
    png = tempname() * ".png"
    @test isfile(plot_reproduction(rep, png)) && filesize(png) > 1000
    # a simulation with failed points (NaN) still plots, with a break instead of a bridge
    gaprun = run_curves(_ -> (xs, [i in 20:25 ? NaN : 2.0 for i in eachindex(xs)]), [h], "gap")
    @test isfile(plot_reproduction(Reproduction(id = "gap", citation = "c", doi = "d", published = gaprun), tempname() * ".png"))
end

# ---- plots -----------------------------------------------------------------------------------

@testset "population plots produce files" begin
    m = onecmt_model()
    s = cohort(8; times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0])
    res = simulate(m, s; replicates = 5, seed = 1)
    dir = mktempdir()
    f1 = spaghetti(res, :C; path = joinpath(dir, "sub", "spaghetti.png"))
    @test f1 isa QSPAtlas.Figure
    @test filesize(joinpath(dir, "sub", "spaghetti.png")) > 1000
    f2 = bands(res, :C; q = (0.05, 0.25, 0.5, 0.75, 0.95), path = joinpath(dir, "bands.png"), logy = true)
    @test f2 isa QSPAtlas.Figure && filesize(joinpath(dir, "bands.png")) > 1000
    @test bands(res, :C; q = (0.1, 0.9), bins = [0.0, 3.0, 12.0], value = :ipred) isa QSPAtlas.Figure       # two q: one band, no line
    @test spaghetti(res, :C; replicate = nothing, value = :ipred) isa QSPAtlas.Figure                     # no path: nothing saved
    obs = [(id = "o$i", time = 1.0 + i, dv = 4.0 / i) for i in 1:5]
    v = vpc(res; endpoint = :C, bins = [0.0, 1.5, 3.0, 6.0, 12.0], observed = obs, q = (0.1, 0.5, 0.9))
    f3 = plot_vpc(v; observed = obs, path = joinpath(dir, "vpc.png"))
    @test f3 isa QSPAtlas.Figure && filesize(joinpath(dir, "vpc.png")) > 1000
    @test plot_vpc(vpc(res; endpoint = :C); logy = true, path = joinpath(dir, "vpc_log.pdf")) isa QSPAtlas.Figure
    @test isfile(joinpath(dir, "vpc_log.pdf"))
    # a VPC with undefined bins and failed values still plots
    badr = simresult([row(k, s, tm, (k == 2 || tm == 3.0) ? NaN : 1.0 + k + tm) for k in 1:12 for s in 1:3 for tm in (1.0, 2.0, 3.0)])
    vnan = @test_logs (:warn, r"non-finite") vpc(badr; endpoint = :C)
    @test plot_vpc(vnan; path = joinpath(dir, "vpc_nan.png")) isa QSPAtlas.Figure
    @test_throws ArgumentError spaghetti(res, :nope)
end

# ---- population safeguards: failure-rate gate, MC precision, tied bootstrap ---------------------------

@testset "failure rate per bin and overall; the gate only marks" begin
    # bin 1: 1 of 5 values failed; bin 2: none; bin 3: 3 of 4
    cols = Dict(1 => [1.0, 2.0, NaN, 4.0, 10.0], 2 => [5.0, 6.0, 7.0, 8.0, 9.0], 3 => [NaN, NaN, 12.0, NaN, missing])
    obs = [row(1, s, tm, ismissing(v) ? NaN : v) for (tm, vs) in cols for (s, v) in enumerate(vs)]
    # the fifth value of bin 3 is NaN too: 4 of 5
    r = simresult(obs)
    p = @test_logs (:warn, r"non-finite") percentiles(r, :C; q = (0.1, 0.5))
    @test [x.failure_rate for x in p if x.q == 0.5] == [0.2, 0.0, 0.8]
    @test [x.n_finite for x in p if x.q == 0.5] == [4, 5, 1]
    @test all(x -> !x.unreliable, p)                                 # default 0.0: report only
    @test failure_rate(p) ≈ 5 / 15                                    # each bin once, not once per quantile
    # the values themselves are as before
    @test at(p; t_mid = 1.0, q = 0.5).value ≈ 3.0 && at(p; t_mid = 3.0, q = 0.5).value == 12.0
    # the gate marks bins above the limit (strictly), does not touch values, warns once
    pg = @test_logs (:warn, r"non-finite") (:warn, r"1 bin\(s\).*unreliable") percentiles(r, :C; q = (0.1, 0.5), max_failure_rate = 0.5)
    @test [x.unreliable for x in pg if x.q == 0.5] == [false, false, true]
    @test [x.value for x in pg] == [x.value for x in p]
    pg2 = @test_logs (:warn, r"non-finite") (:warn, r"2 bin") percentiles(r, :C; q = (0.5,), max_failure_rate = 0.2 - 1e-9)
    @test [x.unreliable for x in pg2] == [true, false, true]
    pg3 = @test_logs (:warn, r"non-finite") percentiles(r, :C; q = (0.5,), max_failure_rate = 0.8)   # 0.8 is not above 0.8
    @test !any(x -> x.unreliable, pg3)
    # an empty bin has no rate and is never marked
    pe = percentiles(r, :C; q = (0.5,), bins = [10.0, 11.0], max_failure_rate = 0.1)
    @test isnan(only(pe).failure_rate) && !only(pe).unreliable && isnan(failure_rate(pe))
    # per replicate rows carry their own rate
    rr = simresult(vcat([row(1, s, 1.0, 1.0) for s in 1:4], [row(2, s, 1.0, s == 1 ? NaN : 1.0) for s in 1:4]))
    pr = @test_logs (:warn, r"non-finite") (:warn, r"unreliable") percentiles(rr, :C; q = (0.5,), by_replicate = true, max_failure_rate = 0.1)
    @test [(x.replicate, x.failure_rate, x.unreliable) for x in pr] == [(1, 0.0, false), (2, 0.25, true)]
    # blq handling does not change the failure rate (it counts non-finite values only)
    pb = @test_logs (:warn, r"non-finite") percentiles(r, :C; q = (0.5,), blq = (lloq = 5.0, handling = :drop))
    @test [x.failure_rate for x in pb] == [0.2, 0.0, 0.8]
    @test_throws ArgumentError percentiles(r, :C; max_failure_rate = -0.1)
    @test_throws ArgumentError percentiles(r, :C; max_failure_rate = 1.5)

    # vpc: per bin over all replicates, and overall
    K = 20
    sim = simresult([row(k, s, tm, (tm == 2.0 && k <= 5) || (tm == 1.0 && k == 1 && s == 1) ? NaN : k + s) for k in 1:K for s in 1:4 for tm in (1.0, 2.0)])
    v = @test_logs (:warn, r"non-finite") vpc(sim; endpoint = :C, q = (0.5,))
    @test [x.failure_rate for x in v.rows] ≈ [1 / 80, 20 / 80]
    @test all(x -> !x.unreliable, v.rows) && v.max_failure_rate == 0.0
    @test v.failure_rate ≈ 21 / 160 && failure_rate(v) == v.failure_rate
    vg = @test_logs (:warn, r"non-finite") (:warn, r"1 bin\(s\).*unreliable") vpc(sim; endpoint = :C, q = (0.5,), max_failure_rate = 0.1)
    @test [x.unreliable for x in vg.rows] == [false, true] && vg.max_failure_rate == 0.1
    @test [x.median for x in vg.rows] == [x.median for x in v.rows]        # values untouched
    @test_throws ArgumentError vpc(sim; endpoint = :C, max_failure_rate = 2.0)
    @test !any(isnan, [vpc(simresult([row(k, 1, 1.0, k) for k in 1:K]); endpoint = :C, q = (0.5,)).failure_rate])
    @test vpc(simresult([row(k, 1, 1.0, k) for k in 1:K]); endpoint = :C, q = (0.5,)).failure_rate == 0.0
    # the old 11-argument VPCResult constructor still works
    old = VPCResult(v.rows, nothing, :C, [0.5], 0.95, K, 0, :dv, nothing, 0, 0)
    @test isnan(old.failure_rate) && old.max_failure_rate == 0.0
end

@testset "mc_check: Monte Carlo precision against a tolerance" begin
    K = 400
    rng = Xoshiro(3)
    sim = simresult([row(k, 1, tm, exp(0.3 * randn(rng))) for k in 1:K for tm in (1.0, 2.0)])
    v = vpc(sim; endpoint = :C, q = (0.1, 0.5, 0.9), seed = 1)
    se = [x.mc_se for x in v.rows]
    @test all(isfinite, se) && length(mc_check(v, 1.0)) == length(v) == 6
    c = mc_check(v, 1.0)
    @test [(x.bin, x.q) for x in c] == [(x.bin, x.q) for x in v.rows] && [x.mc_se for x in c] == se
    @test all(x -> x.limit == 0.2 && x.ok == (x.mc_se < 0.2), c) && all(x -> x.ok, c)
    # the stated fraction of the tolerance: se is about 0.02 here
    tol = maximum(se) / 0.2                                          # limit = max se: strict inequality fails for it
    @test !all(x -> x.ok, mc_check(v, tol)) && count(x -> x.ok, mc_check(v, tol)) == 5
    @test all(x -> x.ok, mc_check(v, 1.001 * tol))
    @test all(x -> x.limit == 0.5 * 2.0, mc_check(v, 2.0; fraction = 0.5))
    # a function gives a tolerance per row, here a fraction of the median
    cf = mc_check(v, r -> 0.5 * r.median; fraction = 0.2)
    @test all(x -> x.ok, cf) && all(isapprox(x.limit, 0.1 * r.median) for (x, r) in zip(cf, v.rows))
    @test !all(x -> x.ok, mc_check(v, 1e-4))
    # fewer than 10 replicates: mc_se is NaN, which is never ok
    few = vpc(simresult([row(k, 1, 1.0, Float64(k)) for k in 1:5]); endpoint = :C, q = (0.5,))
    @test isnan(only(few.rows).mc_se) && !only(mc_check(few, 1e9)).ok
    @test_throws ArgumentError mc_check(v, 0.0)
    @test_throws ArgumentError mc_check(v, -1.0)
    @test_throws ArgumentError mc_check(v, 1.0; fraction = 0.0)
    @test_throws ArgumentError mc_check(v, 1.0; fraction = 1.5)
    @test_throws ArgumentError mc_check(v, r -> "x")
end

@testset "bootstrap SE with heavily tied samples is NaN, not 0" begin
    se(x; seed = 1, kw...) = QSPAtlas.mc_se_median(x, Xoshiro(seed); kw...)
    tied = [zeros(97); ones(3)]                                      # not constant, but every resampled median is 0
    @test length(unique(tied)) == 2
    x = @test_logs (:warn, r"not all equal.*NaN") se(tied)
    @test isnan(x)
    @test @test_logs(se(tied; warn = false)) |> isnan                # silent on request
    # a constant sample really has no error; and the other cases are as before
    @test se(ones(50)) == 0.0
    @test isnan(se(ones(9)))
    @test se([zeros(61); ones(39)]) > 0.05
    # two distinct medians that both occur (balanced 0/1 sample, odd size) are a usable estimate
    bal = [zeros(50); ones(51)]
    xb = @test_logs se(bal)
    @test 0.3 < xb < 0.7
    # the rare one-sided case in the middle: 90 % zeros still gives distinct medians
    @test isfinite(@test_logs se([zeros(55); ones(45)]))
    # the same through vpc: replicate quantiles are 0 in 97 of 100 replicates; one warning for the cell, mc_se NaN
    sim = simresult([row(k, 1, 1.0, k <= 3 ? 1.0 : 0.0) for k in 1:100])
    v = @test_logs (:warn, r"mc_se is NaN in 1 \(bin, q\) cell") vpc(sim; endpoint = :C, q = (0.5,))
    @test isnan(only(v.rows).mc_se) && only(v.rows).median == 0.0
    @test !only(mc_check(v, 1.0)).ok                                 # an unknown precision does not pass
    # all replicate quantiles equal: SE 0 and no warning
    same = simresult([row(k, 1, 1.0, 5.0) for k in 1:100])
    @test only(@test_logs(vpc(same; endpoint = :C, q = (0.5,))).rows).mc_se == 0.0
end


# ---- censored VPC, prediction-corrected VPC, strata, weighting ---------------------------------

# a hand-built result whose subject rows carry covariates (for strata)
function simresult_cov(obs, covs::Dict)
    reps = sort(unique(o.replicate for o in obs))
    subs = [(replicate = r, id = id, η = (;), covariates = covs[id], retcode = :Success) for r in reps for id in sort(collect(keys(covs)))]
    return SimResult(obs, subs, (;))
end

@testset "VPC extensions: censored, prediction-corrected, strata, weighting" begin

@testset "censored quantile follows the type 7 rule and flags the censored ones" begin
    rng = Xoshiro(5)
    for n in (1, 2, 7, 20), q in (0.0, 0.1, 0.25, 0.5, 0.9, 1.0)
        x = randn(rng, n)
        v, below = QSPAtlas.censored_quantile(sort(x), 0, q)
        @test !below && v ≈ quantile(x, q)                         # without censoring it is quantile()
    end
    # 5 values, the lowest 2 censored: order statistics 3, 4, 5 are known (5, 6, 7)
    unc = [5.0, 6.0, 7.0]
    cq(p) = QSPAtlas.censored_quantile(unc, 2, p)
    @test cq(0.2)[2] && isnan(cq(0.2)[1])                          # h = 1.8: touches order statistic 1
    @test cq(0.25)[2]                                              # h = 2: order statistic 2 itself
    @test cq(0.4)[2]                                               # h = 2.6: floor is 2, censored
    @test cq(0.5) == (5.0, false)                                  # h = 3
    @test cq(0.75) == (6.0, false)
    @test cq(0.9)[1] ≈ 6.6 && !cq(0.9)[2]
    @test isnan(QSPAtlas.censored_quantile(Float64[], 0, 0.5)[1]) && !QSPAtlas.censored_quantile(Float64[], 0, 0.5)[2]

    # percentiles: five values at one time, two below the limit
    r = simresult([row(1, s, 1.0, v) for (s, v) in enumerate([0.2, 0.5, 5.0, 6.0, 7.0])])
    p = percentiles(r, :C; q = (0.2, 0.5, 0.75), blq = (lloq = 1.0, handling = :censored))
    @test [(isnan(x.value), x.below_lloq) for x in p] == [(true, true), (false, false), (false, false)]
    @test p[2].value == 5.0 && p[3].value == 6.0 && all(x -> x.n == 5 && x.n_blq == 2, p)
    # other handlings are unchanged: no flag field, a number everywhere
    pk = percentiles(r, :C; q = (0.2,), blq = (lloq = 1.0, handling = :censor_to_lloq))
    @test !haskey(pk[1], :below_lloq) && pk[1].value ≈ 1.0
    @test !haskey(percentiles(r, :C; q = (0.2,))[1], :below_lloq)

    # vpc: every replicate has the same 5 values: the q = 0.2 and 0.4 quantiles are always below the limit
    rb = simresult([row(k, s, 1.0, v) for k in 1:12 for (s, v) in enumerate([0.2, 0.5, 5.0, 6.0, 7.0])])
    v = vpc(rb; endpoint = :C, q = (0.2, 0.5), blq = (lloq = 1.0, handling = :censored))
    x2, x5 = v.rows
    @test isnan(x2.median) && x2.below_lloq && x2.lo_below && x2.hi_below && x2.n_cens == 12 && x2.n_rep == 12 && isnan(x2.mc_se)
    @test x5.median == 5.0 && !x5.below_lloq && x5.n_cens == 0 && x5.n_rep == 12 && x5.lo == 5.0 && x5.hi == 5.0
    @test x2.n_blq == 24
    # replicates that differ: 4 of 10 have the q = 0.4 quantile below the limit, 6 give 3.6 + 0.1 k
    mix = simresult([row(k, s, 1.0, k <= 4 ? v : v2 + 0.1k) for k in 1:10 for (s, (v, v2)) in enumerate(zip([0.2, 0.5, 5.0, 6.0, 7.0], [2.0, 3.0, 4.0, 5.0, 6.0]))])
    vm = vpc(mix; endpoint = :C, q = (0.4,), blq = (lloq = 1.0, handling = :censored), pi_level = 0.5)
    xm = only(vm.rows)
    @test xm.n_cens == 4 && xm.n_rep == 10
    @test xm.median ≈ 4.15 && !xm.below_lloq                       # h = 5.5 over [-, -, -, -, 4.1, ..., 4.6]
    @test xm.lo_below && isnan(xm.lo)                               # the 25 % point (h = 3.25) is among the censored ones
    @test !xm.hi_below && xm.hi ≈ 4.375                              # h = 7.75: between 4.3 and 4.4
    @test isfinite(xm.mc_se)
    # more than half the replicates below the limit: the median is below the limit
    mix6 = simresult([row(k, s, 1.0, k <= 6 ? v : v2 + 0.1k) for k in 1:10 for (s, (v, v2)) in enumerate(zip([0.2, 0.5, 5.0, 6.0, 7.0], [2.0, 3.0, 4.0, 5.0, 6.0]))])
    x6 = only(vpc(mix6; endpoint = :C, q = (0.4,), blq = (lloq = 1.0, handling = :censored)).rows)
    @test x6.below_lloq && isnan(x6.median) && isnan(x6.mc_se)
    # existing handlings do not get the new fields
    @test !haskey(only(vpc(rb; endpoint = :C, q = (0.5,), blq = (lloq = 1.0, handling = :keep)).rows), :below_lloq)
    # an explicit blq flag on observed rows (records stored at the limit)
    obs = [(id = "a$i", time = 1.0, dv = i <= 2 ? 1.0 : 5.0, blq = i <= 2) for i in 1:5]
    vo = vpc(rb; endpoint = :C, q = (0.2, 0.5), blq = (lloq = 1.0, handling = :censored), observed = obs)
    @test vo.observed[1].below_lloq && isnan(vo.observed[1].value) && vo.observed[2].value == 5.0 && vo.observed[1].n_blq == 2
    @test only(vo.blq_fraction).n_observed_blq == 2
end

# five values per time and replicate
rb_small() = simresult([row(k, s, 1.0, v) for k in 1:5 for (s, v) in enumerate([0.5, 0.8, 2.0, 3.0, 4.0])])

@testset "BLQ fraction against a closed form" begin
    # log-normal concentrations, mu = 0, sigma = 1: LLOQ = exp(probit(0.1)) has 10 % of the values below it
    lloq = exp(probit(0.1))
    K, N = 200, 500
    rng = Xoshiro(21)
    sim = simresult([row(k, s, 1.0, exp(randn(rng))) for k in 1:K for s in 1:N])
    obs = [(id = "o$i", time = 1.0, dv = exp(randn(rng))) for i in 1:2000]
    v = vpc(sim; endpoint = :C, q = (0.5,), blq = (lloq = lloq, handling = :censored), observed = obs)
    f = only(v.blq_fraction)
    sd = sqrt(0.1 * 0.9 / N)
    @test f.median ≈ 0.1 atol = 4 * 1.2533 * sd / sqrt(K)
    @test f.lo ≈ 0.1 - 1.96 * sd atol = 0.01
    @test f.hi ≈ 0.1 + 1.96 * sd atol = 0.01
    @test f.n == N && f.n_rep == K && f.mc_se > 0
    @test f.observed ≈ 0.1 atol = 4 * sqrt(0.1 * 0.9 / 2000)
    @test f.n_observed == 2000
    @test !only(v.rows).below_lloq
    @test only(v.rows).median ≈ 1.0 rtol = 0.1
    # no blq, no fraction; an exact small case for the observed fraction
    @test vpc(sim; endpoint = :C, q = (0.5,)).blq_fraction === nothing
    small = vpc(rb_small(), endpoint = :C, q = (0.5,), blq = (lloq = 1.0,), observed = [(id = string(i), time = 1.0, dv = d) for (i, d) in enumerate([0.5, 0.8, 2.0, 3.0, 4.0])])
    fs = only(small.blq_fraction)
    @test fs.observed == 0.4 && fs.n_observed == 5 && fs.n_observed_blq == 2 && fs.median == 0.4 && fs.lo == 0.4 && fs.hi == 0.4
    # the model: ipred of the log-normal clearance model, fraction below the limit known at every time
    m = onecmt_model()
    Nm, Km = 100, 60
    ηL = ω * probit(0.7)                                           # 30 % of the subjects have ipred below L at t = 4
    L = ipred_at(4.0, ηL)
    res = simulate(m, cohort(Nm); replicates = Km, seed = 31, residual = false)
    vm = vpc(res; endpoint = :C, value = :ipred, q = (0.5,), blq = (lloq = L,))
    for fr in vm.blq_fraction
        ηstar = log(-log(L * V0 / DOSE) * V0 / (θCL * fr.t_mid))      # ipred < L  <=>  η > ηstar
        p = 1 - Φ(ηstar / ω)
        @test fr.median ≈ p atol = 4 * 1.2533 * sqrt(p * (1 - p) / Nm) / sqrt(Km) + 1 / Nm
    end
    @test only(f for f in vm.blq_fraction if f.t_mid == 4.0).median ≈ 0.3 atol = 0.04
end

@testset "pcVPC: exact correction on hand-built data" begin
    # six subjects, doses 1 and 5, PRED = dose * 2 at time 1: PRED_bin = median(2, 2, 2, 10, 10, 10) = 6;
    # a dose-proportional value dose * g has pcY = 3 g for every subject, whatever its dose
    K = 20
    dose = Dict("s$s" => (s <= 3 ? 1.0 : 5.0) for s in 1:6)
    g(k, s) = k + 0.1 * s
    sim = simresult([row(k, "s$s", 1.0, dose["s$s"] * g(k, s)) for k in 1:K for s in 1:6])
    table = Dict(("s$s", 1.0) => 2.0 * dose["s$s"] for s in 1:6)
    v = vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = table)
    x = only(v.rows)
    @test v.pc === :pc && v.pred_bin == [6.0] && v.n_nopred == 0
    @test x.median ≈ 3 * (10.5 + 0.35)                             # per replicate 3 (k + 0.35), then the median over k
    @test x.lo ≈ 3 * (quantile(1:K, 0.025) + 0.35) && x.hi ≈ 3 * (quantile(1:K, 0.975) + 0.35)
    # without the correction the doses mix: the values are 1 g and 5 g
    @test only(vpc(sim; endpoint = :C, q = (0.5,)).rows).median != x.median
    # the predictions can also be given as a SimResult of the zero-random-effect run
    predres = SimResult([row(1, "s$s", 1.0, 2.0 * dose["s$s"]) for s in 1:6], [(replicate = 1, id = "s1", η = (;), covariates = (;), retcode = :Success)], (;))
    @test only(vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = predres).rows).median ≈ x.median
    # observed data are corrected with the same PRED_bin (their own prediction through `pred`, or in the row)
    obs = [(id = "s$s", time = 1.0, dv = dose["s$s"] * (5.0 + s), pred = 2.0 * dose["s$s"]) for s in 1:6]
    vo = vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = table, observed = obs)
    @test only(vo.observed).value ≈ 3 * (5.0 + 3.5)
    vo2 = vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = table, observed = [(id = o.id, time = o.time, dv = o.dv) for o in obs])
    @test only(vo2.observed).value ≈ only(vo.observed).value
    # a prediction that is not positive cannot be corrected: left out, counted, warned about
    table0 = merge(table, Dict(("s6", 1.0) => 0.0))
    v0 = @test_logs (:warn, r"no positive finite prediction") vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = table0)
    @test v0.n_nopred == K && only(v0.rows).n == 5
    # PRED_bin is per bin and uses the median of the design (each subject and time once)
    sim2 = simresult(vcat(sim.obs, [row(k, "s$s", 2.0, 100.0 * dose["s$s"]) for k in 1:K for s in 1:6]))
    table2 = merge(table, Dict(("s$s", 2.0) => 20.0 * dose["s$s"] for s in 1:6))
    v2 = vpc(sim2; endpoint = :C, q = (0.5,), pc = true, pred = table2)
    @test v2.pred_bin == [6.0, 60.0]
    @test v2.rows[2].median ≈ 100.0 / 20.0 * 60.0 rtol = 1e-12            # value 100 d, PRED 20 d, PRED_bin 60
    # errors
    @test_throws ArgumentError vpc(sim; endpoint = :C, pc = true)                                   # no predictions
    @test_throws ArgumentError vpc(sim; endpoint = :C, pc = true, pred = Dict(("s1", 1.0) => 2.0))    # incomplete predictions
    @test_throws ArgumentError vpc(sim; endpoint = :C, variability_correction = true)               # needs pc
    @test_throws ArgumentError vpc(simresult([row(1, s, 1.0, 1.0) for s in 1:3]); endpoint = :C, pc = true, variability_correction = true,
        pred = Dict((string(s), 1.0) => 1.0 for s in 1:3))                                          # one replicate
    @test_throws ArgumentError vpc(sim; endpoint = :C, pc = true, pred = 3)

    # variability correction: values dose * s * k, so the spread of observation s grows with s
    sim3 = simresult([row(k, "s$s", 1.0, dose["s$s"] * s * k) for k in 1:K for s in 1:6])
    v3 = vpc(sim3; endpoint = :C, q = (0.5,), pc = true, variability_correction = true, pred = table)
    @test v3.pc === :pvc
    pcy(k, s) = dose["s$s"] * s * k * 6.0 / (2.0 * dose["s$s"])                                     # 3 s k
    sdij(s) = std([pcy(k, s) for k in 1:K])
    sdbin = median([sdij(s) for s in 1:6])
    pvc(k, s) = 6.0 + (pcy(k, s) - 6.0) * sdbin / sdij(s)
    @test only(v3.rows).median ≈ median([median([pvc(k, s) for s in 1:6]) for k in 1:K])
    # with it the spread over subjects is the same as the spread over replicates would imply: not the pc one
    @test only(v3.rows).median != only(vpc(sim3; endpoint = :C, q = (0.5,), pc = true, pred = table).rows).median
    # an observed record outside the simulated design cannot be variability-corrected
    @test_throws ArgumentError vpc(sim3; endpoint = :C, q = (0.5,), pc = true, variability_correction = true, pred = table,
        observed = [(id = "s1", time = 1.0, dv = 1.0, pred = 2.0), (id = "zz", time = 1.0, dv = 1.0, pred = 2.0)])
end

@testset "pcVPC and strata on a dose-proportional model" begin
    # linear PK, ipred = dose / V exp(-CL e^η t / V): ipred / PRED = exp(-CL (e^η - 1) t / V) for every dose
    m = onecmt_model()
    N1, K = 60, 60
    doses = [100.0, 500.0]
    subs = [Subject("d$(Int(d))_$i"; regimen = Regimen(d; cmt = :central), covariates = (dose = d,), obstimes = [1.0, 4.0, 12.0]) for d in doses for i in 1:N1]
    res = simulate(m, subs; replicates = K, seed = 41, residual = false)
    zrun = simulate(m, subs; η = :zero, residual = false)
    # the corrected quantile in units of PRED_bin: closed form, and the finite-N slack of the sample quantile (+-1/N in probability, as above)
    ratio(tm, qq) = exp(-θCL * (exp(ω * probit(1 - qq)) - 1) * tm / V0)
    slack(x, pb, tm, qq, N) = 4 * x.mc_se / pb + abs(ratio(tm, clamp(qq - 1 / N, 0, 1)) - ratio(tm, clamp(qq + 1 / N, 0, 1)))
    # (a) unstratified pcVPC of the mixed-dose cohort
    pooled = vpc(res; endpoint = :C, value = :ipred, pc = true, pred = zrun)
    zre(s, tm) = only(o.ipred for o in zero_re_prediction(m, s).obs if o.time == tm)       # PRED of one observation
    for (b, tm) in enumerate((1.0, 4.0, 12.0))
        @test pooled.pred_bin[b] ≈ median([zre(s, tm) for s in subs]) rtol = 1e-8
    end
    for x in pooled.rows
        pb = pooled.pred_bin[x.bin]
        @test abs(x.median / pb - ratio(x.t_mid, x.q)) <= slack(x, pb, x.t_mid, x.q, 2 * N1)
    end
    raw = vpc(res; endpoint = :C, value = :ipred)
    # without the correction the doses mix: the 90 % over 10 % quantile ratio is much wider
    spread(v, tm) = at(v.rows; t_mid = tm, q = 0.9).median / at(v.rows; t_mid = tm, q = 0.1).median
    @test spread(raw, 4.0) > 2 * spread(pooled, 4.0)
    # (b) stratified by dose: the strata collapse after the correction, differ by the dose ratio before it
    by_dose = s -> s.covariates.dose
    sraw = vpc(res; endpoint = :C, value = :ipred, strata = by_dose)
    spc = vpc(res; endpoint = :C, value = :ipred, strata = by_dose, pc = true, pred = zrun)
    @test sraw isa StratifiedVPC && spc isa StratifiedVPC && collect(keys(spc)) == doses
    @test all(x -> x.n == N1, spc[100.0].rows) && all(x -> x.n == N1, spc[500.0].rows)
    for (b, tm) in enumerate((1.0, 4.0, 12.0)), qq in (0.1, 0.5, 0.9)
        r1, r5 = at(sraw[100.0].rows; t_mid = tm, q = qq), at(sraw[500.0].rows; t_mid = tm, q = qq)
        @test r5.median / r1.median ≈ 5.0 rtol = 0.1                    # dose-proportional, before the correction
        p1, p5 = at(spc[100.0].rows; t_mid = tm, q = qq), at(spc[500.0].rows; t_mid = tm, q = qq)
        pb1, pb5 = spc[100.0].pred_bin[b], spc[500.0].pred_bin[b]
        @test pb5 / pb1 ≈ 5.0 rtol = 1e-8
        @test abs(p1.median / pb1 - ratio(tm, qq)) <= slack(p1, pb1, tm, qq, N1)
        @test abs(p5.median / pb5 - ratio(tm, qq)) <= slack(p5, pb5, tm, qq, N1)
        @test abs(p5.median / pb5 - p1.median / pb1) <= slack(p1, pb1, tm, qq, N1) + slack(p5, pb5, tm, qq, N1)   # collapsed
    end
    # the model method computes the predictions itself and gives the same result; observed data at a time that is
    # not a cohort observation time get their prediction from the model
    obs = [(id = s.id, time = 2.5, dv = 3.0) for s in subs[1:5]]
    kw = (endpoint = :C, replicates = 10, seed = 3, value = :ipred, bins = [0.0, 3.0, 6.0, 12.0], observed = obs, q = (0.5,))
    vo = vpc(m, subs; kw..., pc = true)
    @test vo.n_unbinned == 0 && vo.observed[1].n == 5 && vo.pc === :pc && vo.n_nopred == 0
    pred25 = [only(o.ipred for o in simulate(m, Subject(s.id, s.regimen, s.covariates, s.covariate_times, s.covariate_values, [2.5], s.control); η = :zero, residual = false).obs) for s in subs[1:5]]
    @test vo.observed[1].value ≈ median(3.0 .* vo.pred_bin[1] ./ pred25)
    res10 = simulate(m, subs; replicates = 10, seed = 3)
    kw0 = Base.structdiff(kw, (; replicates = 0))
    vref = vpc(res10; kw0..., pc = true, pred = merge(Dict((o.id, o.time) => o.ipred for o in zrun.obs),
        Dict((s.id, 2.5) => p for (s, p) in zip(subs[1:5], pred25))))
    @test isequal(vo.rows, vref.rows) && isequal(vo.observed, vref.observed)
    # strata and the pcVPC together through the model method; variability correction needs >= 2 replicates
    vs = vpc(m, subs; endpoint = :C, replicates = 10, seed = 3, value = :ipred, q = (0.5,), pc = true, variability_correction = true, strata = by_dose)
    @test all(v -> v.pc === :pvc, vs.results)
    @test vs[500.0].pred_bin[1] / vs[100.0].pred_bin[1] ≈ 5.0 rtol = 1e-8
    @test_throws ArgumentError vpc(m, subs; endpoint = :C, replicates = 1, pc = true, variability_correction = true)
    @test_throws ArgumentError vpc(m, subs; endpoint = :C, replicates = 5, variability_correction = true)
end

@testset "stratification: one result per stratum, counts and observed data" begin
    covs = Dict("a1" => (grp = "A",), "a2" => (grp = "A",), "b1" => (grp = "B",), "b2" => (grp = "B",), "b3" => (grp = "B",))
    K = 12
    sim = simresult_cov([row(k, id, tm, k + 10 * (id[1] == 'b') + tm) for k in 1:K for id in keys(covs) for tm in (1.0, 2.0)], covs)
    sv = vpc(sim; endpoint = :C, q = (0.5,), strata = s -> s.covariates.grp)
    @test sv isa StratifiedVPC && length(sv) == 2 && collect(keys(sv)) == ["A", "B"]
    @test occursin("StratifiedVPC", sprint(show, sv))
    @test all(x -> x.n == 2, sv["A"].rows) && all(x -> x.n == 3, sv["B"].rows)
    @test sv["A"].stratum == "A" && sv["B"].stratum == "B" && sv[2] === sv["B"]
    @test [k for (k, _) in sv] == ["A", "B"] && [v.stratum for (_, v) in sv] == ["A", "B"]
    @test all(x -> x.n_rep == K, sv["A"].rows)
    # each stratum is the unstratified VPC of its own subjects
    subA = SimResult(filter(o -> o.id in ("a1", "a2"), sim.obs), filter(s -> s.id in ("a1", "a2"), sim.subjects), sim.settings)
    @test isequal(vpc(subA; endpoint = :C, q = (0.5,)).rows, sv["A"].rows)
    @test_throws KeyError sv["C"]
    # the dictionary form, same result
    sd = vpc(sim; endpoint = :C, q = (0.5,), strata = Dict(id => c.grp for (id, c) in covs))
    @test isequal(sd["B"].rows, sv["B"].rows)
    # observed data are split by the stratum of their subject
    obs = vcat([(id = "a1", time = 1.0, dv = 2.0), (id = "a2", time = 1.0, dv = 4.0)], [(id = "b$i", time = 1.0, dv = 10.0 + i) for i in 1:3])
    svo = vpc(sim; endpoint = :C, q = (0.5,), strata = s -> s.covariates.grp, observed = obs)
    @test at(svo["A"].observed; t_mid = 1.0).value == 3.0 && at(svo["A"].observed; t_mid = 1.0).n == 2
    @test at(svo["B"].observed; t_mid = 1.0).value == 12.0 && at(svo["B"].observed; t_mid = 1.0).n == 3
    @test at(svo["A"].observed; t_mid = 2.0).n == 0
    @test_throws ArgumentError vpc(sim; endpoint = :C, strata = s -> s.covariates.grp, observed = [(id = "zz", time = 1.0, dv = 1.0)])
    # no strata: unchanged type
    @test vpc(sim; endpoint = :C, q = (0.5,)) isa VPCResult && vpc(sim; endpoint = :C, q = (0.5,)).stratum === nothing
    # the model method passes the Subject: split by a covariate, the simulation is shared
    m = onecmt_model()
    sub = vcat([Subject("lo$i"; regimen = Regimen(DOSE; cmt = :central), covariates = (arm = :low,), obstimes = [1.0, 4.0]) for i in 1:3],
        [Subject("hi$i"; regimen = Regimen(2 * DOSE; cmt = :central), covariates = (arm = :high,), obstimes = [1.0, 4.0]) for i in 1:5])
    sm = vpc(m, sub; endpoint = :C, replicates = 10, seed = 2, strata = s -> s.covariates.arm)
    @test collect(keys(sm)) == [:high, :low] && all(x -> x.n == 5, sm[:high].rows) && all(x -> x.n == 3, sm[:low].rows)
    whole = vpc(m, sub; endpoint = :C, replicates = 10, seed = 2)
    @test all(x -> x.n == 8, whole.rows)
    # plots, per stratum, with and without the BLQ panel
    dir = mktempdir()
    @test plot_vpc(sv; path = joinpath(dir, "strata.png")) isa QSPAtlas.Figure && filesize(joinpath(dir, "strata.png")) > 1000
    svb = vpc(sim; endpoint = :C, q = (0.1, 0.5, 0.9), strata = s -> s.covariates.grp, blq = (lloq = 3.0, handling = :censored), observed = obs)
    @test plot_vpc(svb; blq_fraction = true, ncols = 1, path = joinpath(dir, "strata_blq.png")) isa QSPAtlas.Figure
    @test_throws ArgumentError plot_vpc(sv; blq_fraction = true)
    @test_throws ArgumentError plot_vpc(sv; ncols = 0)
end

@testset "BLQ panel of plot_vpc" begin
    rb = simresult([row(k, s, tm, v * tm) for k in 1:12 for (s, v) in enumerate([0.2, 0.5, 5.0, 6.0, 7.0]) for tm in (1.0, 2.0, 3.0)])
    obs = [(id = "o$i", time = 1.0 + i % 3, dv = 0.3 * i) for i in 1:15]
    v = vpc(rb; endpoint = :C, blq = (lloq = 1.0, handling = :censored), observed = obs)
    dir = mktempdir()
    @test plot_vpc(v; blq_fraction = true, path = joinpath(dir, "blq.png")) isa QSPAtlas.Figure && filesize(joinpath(dir, "blq.png")) > 1000
    @test plot_vpc(v; blq_fraction = true, logy = true) isa QSPAtlas.Figure
    @test plot_vpc(vpc(rb; endpoint = :C, blq = (lloq = 1.0,)); blq_fraction = true) isa QSPAtlas.Figure     # no observed data
    @test_throws ArgumentError plot_vpc(vpc(rb; endpoint = :C); blq_fraction = true)
    # a prediction-corrected result plots too
    tbl = Dict((string(s), tm) => 1.0 + s for s in 1:5 for tm in (1.0, 2.0, 3.0))
    @test plot_vpc(vpc(rb; endpoint = :C, pc = true, pred = tbl)) isa QSPAtlas.Figure
end

@testset "weighting: observation (pooled) or subject" begin
    # subject A has four observations in the bin (1, 2, 3, 4 at times 1..4), B and C one each (10, 10)
    mk(k) = vcat([row(k, "A", tm, tm) for tm in 1.0:4.0], [row(k, "B", 2.0, 10.0), row(k, "C", 3.0, 10.0)])
    sim = simresult(vcat([mk(k) for k in 1:10]...))
    edges = [0.0, 10.0]
    byobs = only(percentiles(sim, :C; q = (0.5,), bins = edges))
    @test byobs.n == 60 && byobs.value == 3.5                        # pooled: 1, 2, 3, 4, 10, 10 per replicate
    for (sv, a) in ((:mean, 2.5), (:first, 1.0), (:last, 4.0))
        x = only(percentiles(sim, :C; q = (0.5,), bins = edges, weighting = :subject, subject_value = sv))
        @test x.n == 30 && x.n_finite == 60 && x.n_failed == 0       # one value per subject and replicate, n_finite counts observations
        @test x.value == 10.0                                        # subject values (a, 10, 10): median 10
        @test only(percentiles(sim, :C; q = (0.01,), bins = edges, weighting = :subject, subject_value = sv)).value == a   # the subject's own value
    end
    # the smallest subject value is each subject's own, per replicate
    r1 = only(percentiles(sim, :C; q = (0.0,), bins = edges, weighting = :subject, subject_value = :last, by_replicate = true)[1:1])
    @test r1.value == 4.0 && r1.n == 3
    @test only(percentiles(sim, :C; q = (0.0,), bins = edges, weighting = :subject, subject_value = :first, by_replicate = true)[1:1]).value == 1.0
    @test only(percentiles(sim, :C; q = (0.0,), bins = edges, weighting = :subject, subject_value = :mean, by_replicate = true)[1:1]).value == 2.5
    # vpc: the same per replicate, so the spread is zero and the medians are those of the subject values
    vo = only(vpc(sim; endpoint = :C, q = (0.5,), bins = edges).rows)
    vs = only(vpc(sim; endpoint = :C, q = (0.5,), bins = edges, weighting = :subject).rows)
    @test vo.median == 3.5 && vo.n == 6 && vs.median == 10.0 && vs.n == 3
    @test vpc(sim; endpoint = :C, bins = edges, weighting = :subject).weighting === :subject && vpc(sim; endpoint = :C).weighting === :observation
    # observed data follow the same weighting
    obs = vcat([(id = "A", time = tm, dv = tm) for tm in 1.0:4.0], [(id = "B", time = 2.0, dv = 10.0), (id = "C", time = 3.0, dv = 10.0)])
    @test only(vpc(sim; endpoint = :C, q = (0.5,), bins = edges, observed = obs).observed).value == 3.5
    @test only(vpc(sim; endpoint = :C, q = (0.5,), bins = edges, observed = obs, weighting = :subject).observed).value == 10.0
    # one observation per subject and bin: no difference
    sim1 = simresult([row(k, s, 1.0, k + s) for k in 1:8 for s in 1:5])
    @test isequal(vpc(sim1; endpoint = :C).rows, vpc(sim1; endpoint = :C, weighting = :subject).rows)
    # the BLQ fraction also counts subjects with weighting = :subject: A has 1 of 4 values below 1.5, so its mean is censored
    vb = vpc(sim; endpoint = :C, q = (0.5,), bins = edges, blq = (lloq = 1.5, handling = :keep))
    @test only(vb.blq_fraction).median ≈ 1 / 6
    vbs = vpc(sim; endpoint = :C, q = (0.5,), bins = edges, blq = (lloq = 1.5, handling = :keep), weighting = :subject)
    @test only(vbs.blq_fraction).median ≈ 1 / 3
    vbf = vpc(sim; endpoint = :C, q = (0.5,), bins = edges, blq = (lloq = 1.5, handling = :keep), weighting = :subject, subject_value = :first)
    @test only(vbf.blq_fraction).median ≈ 1 / 3
    @test_throws ArgumentError percentiles(sim, :C; weighting = :nope)
    @test_throws ArgumentError vpc(sim; endpoint = :C, subject_value = :median)
end

end

@testset "VPC review fixes: records, censoring limits, subject means, typical predictions" begin

@testset "pcVPC keys predictions by record: a pre-dose and a post-dose record at one time" begin
    # per replicate k and subject s: a pre-dose record s k (PRED 1) and a post-dose record 10 s k (PRED 10), both at time 1
    K = 10
    rec(k, s) = [row(k, "s$s", 1.0, s * k), row(k, "s$s", 1.0, 10.0 * s * k)]
    sim = simresult(vcat([rec(k, s) for k in 1:K for s in 1:4]...))
    pred3 = Dict{Tuple{String, Float64, Int}, Float64}()
    for s in 1:4
        pred3[("s$s", 1.0, 1)] = 1.0
        pred3[("s$s", 1.0, 2)] = 10.0
    end
    v = vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = pred3)
    # PRED_bin is the median over the 8 records (1, 10 four times each): 5.5; every record becomes 5.5 s k
    @test v.pred_bin == [5.5] && v.n_nopred == 0
    @test only(v.rows).median ≈ median([median([5.5 * s * k for s in 1:4 for _ in 1:2]) for k in 1:K])
    @test only(v.rows).n == 8
    # the same from a SimResult of the typical predictions, which has both records
    predres = SimResult(vcat([[row(1, "s$s", 1.0, 1.0), row(1, "s$s", 1.0, 10.0)] for s in 1:4]...),
        [(replicate = 1, id = "s$s", η = (;), covariates = (;), retcode = :Success) for s in 1:4], (;))
    @test isequal(vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = predres).rows, v.rows)
    # a zero pre-dose prediction leaves out the pre-dose records only, never the post-dose ones
    pz = merge(pred3, Dict(("s$s", 1.0, 1) => 0.0 for s in 1:4))
    vz = @test_logs (:warn, r"no positive finite prediction") vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = pz)
    @test vz.n_nopred == 4K && only(vz.rows).n == 4 && vz.pred_bin == [10.0]
    @test only(vz.rows).median ≈ median([median([10.0 * s * k for s in 1:4]) for k in 1:K])
    # observed records are numbered the same way, `order = :pre` first whatever the row order
    obs = [(id = "s1", time = 1.0, dv = 1.0, order = :post), (id = "s1", time = 1.0, dv = 1.0, order = :pre)]
    vo = vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = pred3, observed = obs)
    @test only(vo.observed).value ≈ (5.5 / 1.0 + 5.5 / 10.0) / 2
    @test QSPAtlas.obs_record_keys(obs, :C) == [("s1", 1.0, 2), ("s1", 1.0, 1)]
    # variability correction: sd_ij per record (the post-dose records have twice the spread)
    sim3 = simresult(vcat([[row(k, "s$s", 1.0, s * k), row(k, "s$s", 1.0, 20.0 * s * k)] for k in 1:K for s in 1:4]...))
    v3 = vpc(sim3; endpoint = :C, q = (0.5,), pc = true, variability_correction = true, pred = pred3)
    pcy(k, s, j) = j == 1 ? 5.5 * s * k : 11.0 * s * k
    sdij(s, j) = std([pcy(k, s, j) for k in 1:K])
    sdbin = median([sdij(s, j) for s in 1:4 for j in 1:2])
    pvc(k, s, j) = 5.5 + (pcy(k, s, j) - 5.5) * sdbin / sdij(s, j)
    @test only(v3.rows).median ≈ median([median([pvc(k, s, j) for s in 1:4 for j in 1:2]) for k in 1:K])
    # a dictionary keyed by (id, time) applies to every record at that time
    v2 = vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = Dict(("s$s", 1.0) => 2.0 for s in 1:4))
    @test v2.pred_bin == [2.0]

    # the model method: a second dose at 6, sampled just before and just after it
    m = onecmt_model()
    k_el = θCL / V0
    p1 = DOSE / V0 * exp(-k_el * 6.0)                  # typical pre-dose prediction
    p2 = p1 + DOSE / V0                                # typical post-dose prediction
    reg = Regimen(DOSE; cmt = :central, ii = 6.0, addl = 1)
    subs = [Subject("p$i"; regimen = reg, obstimes = [6.0, 6.0, 10.0], obsorder = [:pre, :post, :post]) for i in 1:10]
    push!(subs, Subject("x"; regimen = reg, obstimes = [10.0]))   # its observed records at 6 are not in the design
    obsm = [(id = "p1", time = 6.0, dv = 1.0, order = :post), (id = "p1", time = 6.0, dv = 1.0, order = :pre),
        (id = "x", time = 6.0, dv = 1.0, order = :pre), (id = "x", time = 6.0, dv = 1.0, order = :post)]
    vm = vpc(m, subs; endpoint = :C, value = :ipred, q = (0.5,), replicates = 4, seed = 7, pc = true, observed = obsm)
    @test vm.pred_bin[1] ≈ (p1 + p2) / 2 rtol = 1e-7
    @test vm.observed[1].value ≈ (vm.pred_bin[1] / p1 + vm.pred_bin[1] / p2) / 2 rtol = 1e-7
    # by hand from the simulation: each record corrected by its own prediction
    res = simulate(m, subs; replicates = 4, seed = 7)
    byrep = map(1:4) do k
        rows6 = [o for o in res.obs if o.replicate == k && o.time == 6.0]
        vals = Float64[]
        for id in unique(o.id for o in rows6)
            r2 = [o for o in rows6 if o.id == id]
            push!(vals, r2[1].ipred * vm.pred_bin[1] / p1, r2[2].ipred * vm.pred_bin[1] / p2)
        end
        median(vals)
    end
    @test vm.rows[1].median ≈ median(byrep) rtol = 1e-7
end

@testset "censored quantiles with a limit per value" begin
    rng = Xoshiro(9)
    for n in (1, 3, 8), nc in (0, 1, 3), q in (0.0, 0.2, 0.5, 0.85, 1.0)
        u = sort(1.0 .+ rand(rng, n))
        a = QSPAtlas.censored_quantile(u, nc, q)
        b = QSPAtlas.censored_quantile_limits(u, fill(1.0, nc), q)
        @test isequal(a[1], b[1]) && a[2] == b[2] && (b[2] ? b[3] == 1.0 : isnan(b[3]))    # a common limit: as before
    end
    cl(u, l, p) = QSPAtlas.censored_quantile_limits(u, l, p)
    # quantified 2, 3, 4 and one value below 3.5: order statistics 1 to 3 depend on it, 4 does not
    @test isequal(cl([2.0, 3.0, 4.0], [3.5], 1.0), (4.0, false, NaN))
    x = cl([2.0, 3.0, 4.0], [3.5], 2 / 3)             # order statistic 3: 3 or 3.5
    @test isnan(x[1]) && x[2] && x[3] == 3.5
    x = cl([2.0, 3.0, 4.0], [5.0], 1.0)                # a limit above every quantified value: nothing is known
    @test isnan(x[1]) && x[2] && x[3] == 5.0
    # ties: with 5, 5 and one value below 6, order statistic 2 is 5 wherever it lies, 3 is not known
    @test isequal(cl([5.0, 5.0], [6.0], 0.5), (5.0, false, NaN))
    @test cl([5.0, 5.0], [6.0], 1.0)[2]
    @test isnan(cl(Float64[], Float64[], 0.5)[1])

    # pcVPC with censoring: A (PRED 1) below LLOQ = 1, B (PRED 10) at 3, C (PRED 10) at 20. PRED_bin = 10:
    # A's limit becomes 10, B's corrected value 3. Ranking A below B would give the median 3.
    K = 12
    sim = simresult([row(k, id, 1.0, v) for k in 1:K for (id, v) in (("A", 0.5), ("B", 3.0), ("C", 20.0))])
    pr = Dict(("A", 1.0) => 1.0, ("B", 1.0) => 10.0, ("C", 1.0) => 10.0)
    obs = [(id = "A", time = 1.0, dv = 0.5), (id = "B", time = 1.0, dv = 3.0), (id = "C", time = 1.0, dv = 20.0)]
    v = vpc(sim; endpoint = :C, q = (0.5, 1.0), pc = true, pred = pr, blq = (lloq = 1.0, handling = :censored), observed = obs)
    @test v.pred_bin == [10.0]
    x5, x1 = v.rows
    @test x5.below_lloq && isnan(x5.median) && x5.n_cens == K && x5.n_blq == K
    @test x1.median == 20.0 && !x1.below_lloq
    @test v.observed[1].below_lloq && isnan(v.observed[1].value) && v.observed[1].limit == 10.0
    @test v.observed[2].value == 20.0 && isnan(v.observed[2].limit)
    # without the correction the common limit gives the usual answer
    @test only(vpc(sim; endpoint = :C, q = (0.5,), blq = (lloq = 1.0, handling = :censored)).rows).median == 3.0
    # :censor_to_lloq puts A at its own corrected limit, 10
    @test only(vpc(sim; endpoint = :C, q = (0.5,), pc = true, pred = pr, blq = (lloq = 1.0, handling = :censor_to_lloq)).rows).median == 10.0
    # replicates whose quantile is below a limit are summarised with their own bounds: in half the
    # replicates A is quantified (2, corrected 20), in the others it is below its limit 10
    simh = simresult([row(k, id, 1.0, id == "A" ? (k <= 6 ? 0.5 : 2.0) : v) for k in 1:K for (id, v) in (("A", 0.5), ("B", 3.0), ("C", 20.0))])
    vh = vpc(simh; endpoint = :C, q = (0.5,), pc = true, pred = pr, blq = (lloq = 1.0, handling = :censored))
    # quantified replicates: median of (20, 3, 20) = 20; the others are below 10 < 20, so the median of
    # (6 x below 10, 6 x 20) has its lower order statistic among the censored ones
    @test only(vh.rows).below_lloq && only(vh.rows).n_cens == 6 && only(vh.rows).hi == 20.0
    # the plot of a prediction-corrected result draws no single LLOQ line; the uncorrected one does
    hl(vv) = begin
        fig = plot_vpc(vv)
        axs = [c for c in fig.content if c isa QSPAtlas.MakieAxis]
        any(p -> occursin("hlines", lowercase(string(typeof(p)))), only(axs).scene.plots)
    end
    @test !hl(v)
    @test hl(vpc(sim; endpoint = :C, q = (0.5,), blq = (lloq = 1.0, handling = :censored)))
end

@testset "subject means that include a censored value are censored" begin
    # A: 0.5 (below 1) and 10, B: 3 and 3, C: 4 and 4, all in one bin
    mk(a) = simresult([row(1, "A", 1.0, a), row(1, "A", 2.0, 10.0), row(1, "B", 1.0, 3.0), row(1, "B", 2.0, 3.0),
        row(1, "C", 1.0, 4.0), row(1, "C", 2.0, 4.0)])
    edges = [0.0, 10.0]
    kw = (q = (0.5, 1.0), bins = edges, weighting = :subject)
    p = percentiles(mk(0.5), :C; kw..., blq = (lloq = 1.0, handling = :censored))
    # A's mean is censored with limit (1 + 10) / 2 = 5.5 above B and C: neither the median nor the maximum is known
    @test all(x -> x.below_lloq && isnan(x.value), p) && p[1].limit == 4.0 && p[2].limit == 5.5
    @test all(x -> x.n == 3 && x.n_blq == 1, p)
    # the placeholder of the censored value does not matter
    @test isequal(p, percentiles(mk(0.9), :C; kw..., blq = (lloq = 1.0, handling = :censored)))
    # :drop leaves A out of the bin; :keep uses its number and counts it
    pd = percentiles(mk(0.5), :C; kw..., blq = (lloq = 1.0, handling = :drop))
    @test pd[1].value == 3.5 && pd[1].n == 2 && pd[1].n_blq == 1
    pk = percentiles(mk(0.5), :C; kw..., blq = (lloq = 1.0, handling = :keep))
    @test pk[2].value == 5.25 && pk[1].n_blq == 1
    # :first picks A's censored value, whose limit is the LLOQ: the median is B's
    pf = percentiles(mk(0.5), :C; kw..., subject_value = :first, blq = (lloq = 1.0, handling = :censored))
    @test pf[1].value == 3.0 && !pf[1].below_lloq
end

end

end
