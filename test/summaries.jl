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

end
