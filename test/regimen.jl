# Dosing regimens: expansion, F and lag, units, and solver checks against closed forms.
# Own module, so helper names (onecmt, TOL, ...) cannot clash with other test files.
module RegimenTests

using QSPAtlas, Test
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicIndexingInterface: getp
using SciMLBase: EnsembleProblem, EnsembleThreads, CallbackSet

const TOL = (abstol = 1e-12, reltol = 1e-12)
close(a, b; rtol = 1e-7, atol = 1e-10) = isapprox(a, b; rtol, atol)

# One compartment with an infusion input: dA/dt = -k A + Rin (k = 0.3, Rin = 0 unless set).
function onecmt(; kval = 0.3, rin = 0.0)
    @parameters k = kval Rin = rin
    @variables A(t) = 0.0
    @mtkcompile sys = System([D(A) ~ -k * A + Rin], t)
    return sys
end

# Oral: depot -> central, infusion into central.
function oral(; kaval = 1.1, kval = 0.3)
    @parameters ka = kaval k = kval Rin = 0.0
    @variables Ad(t) = 0.0 Ac(t) = 0.0
    @mtkcompile sys = System([D(Ad) ~ -ka * Ad, D(Ac) ~ ka * Ad - k * Ac + Rin], t)
    return sys
end

# Run a regimen on `sys`; returns the solution evaluated on `ts`.
function run_regimen(sys, routes, regimen, tspan; control = DoseControl(), alg = Tsit5(), ts = nothing, extra = [], tol = TOL)
    ev = resolve(regimen, routes, control)
    cb, tstops = dosing_callback(sys, ev, routes)
    prob = ODEProblem(sys, extra, tspan; callback = cb, tstops = tstops)
    sol = solve(prob, alg; tol...)
    return sol, prob
end

# Exact amount for a one-compartment infusion `rate` on [t1, t2] at time x.
function infusion_exact(x, t1, t2, rate, k)
    x <= t1 && return 0.0
    x <= t2 && return rate / k * (1 - exp(-k * (x - t1)))
    return rate / k * (1 - exp(-k * (t2 - t1))) * exp(-k * (x - t2))
end

@testset "Regimen expansion" begin
    r = Regimen(100.0; cmt = :depot, time = 2.0, ii = 12.0, addl = 2)
    @test [d.time for d in doses(r)] == [2.0, 14.0, 26.0]
    @test all(d -> d.amt == 100.0 && d.cmt === :depot && d.rate == 0.0 && d.duration == 0.0, doses(r))
    @test length(doses(Regimen(5.0; cmt = :a))) == 1 && doses(Regimen(5.0; cmt = :a))[1].time == 0.0

    # hand-written list
    hand = [Dose(0.0, 10.0; cmt = :a), Dose(8.0, 10.0; cmt = :a), Dose(16.0, 10.0; cmt = :a)]
    r2 = Regimen(10.0; cmt = :a, ii = 8.0, addl = 2)
    @test [(d.time, d.amt, d.cmt) for d in doses(r2)] == [(d.time, d.amt, d.cmt) for d in hand]

    # concatenation sorts by time and keeps the order of equal-time doses
    a = Regimen(1.0; cmt = :a, time = 5.0)
    b = Regimen(2.0; cmt = :b, time = 0.0, ii = 5.0, addl = 1)   # times 0, 5
    c = vcat(a, b)
    @test [(d.time, d.cmt) for d in doses(c)] == [(0.0, :b), (5.0, :a), (5.0, :b)]
    @test [d.amt for d in doses(Regimen(a, b))] == [2.0, 1.0, 2.0]
    @test doses(Regimen(b, a)) |> x -> [(d.time, d.cmt) for d in x] == [(0.0, :b), (5.0, :b), (5.0, :a)]
    @test isempty(doses(Regimen()))

    # infusions: rate fixes duration, duration fixes rate
    ir = Dose(1.0, 100.0; cmt = :central, rate = 20.0)
    @test ir.duration == 5.0 && ir.fixed === :rate
    id = Dose(1.0, 100.0; cmt = :central, duration = 4.0)
    @test id.rate == 25.0 && id.fixed === :duration
    both = Dose(0.0, 100.0; cmt = :central, rate = 20.0, duration = 5.0)
    @test both.rate == 20.0 && both.duration == 5.0
    ri = doses(Regimen(100.0; cmt = :central, ii = 24.0, addl = 1, duration = 2.0))
    @test all(d -> d.rate == 50.0 && d.duration == 2.0, ri)
end

@testset "Dosing validation" begin
    @test_throws ArgumentError Dose(0.0, -1.0; cmt = :a)
    @test_throws ArgumentError Dose(0.0, NaN; cmt = :a)
    @test_throws ArgumentError Dose(NaN, 1.0; cmt = :a)
    @test_throws ArgumentError Dose(0.0, 1.0; cmt = :a, rate = -1.0)
    @test_throws ArgumentError Dose(0.0, 1.0; cmt = :a, duration = -1.0)
    @test_throws ArgumentError Dose(0.0, 100.0; cmt = :a, rate = 20.0, duration = 6.0)   # amt/rate = 5
    @test_throws ArgumentError Regimen(1.0; cmt = :a, ii = 0.0, addl = 2)
    @test_throws ArgumentError Regimen(1.0; cmt = :a, ii = 1.0, addl = -1)
    @test_throws ArgumentError Regimen(1.0; cmt = :a, ii = -1.0, addl = 1)
    @test_throws ArgumentError DoseControl(bioav = Dict(:a => -0.1))
    @test_throws ArgumentError DoseControl(lag = Dict(:a => -1.0))
    sys = onecmt()
    @test_throws ArgumentError Route(:central, sys.A; scale = 0.0)
    r = Route(:central, sys.A; rate = sys.Rin)
    @test_throws ArgumentError resolve(Regimen(1.0; cmt = :nowhere), [r])          # no route
    @test_throws ArgumentError resolve(Regimen(1.0; cmt = :central), [r, r])      # duplicate routes
    @test_throws ArgumentError resolve(Regimen(1.0; cmt = :central, rate = 1.0), [Route(:central, sys.A)])   # no rate parameter
    @test resolve(Regimen(1.0; cmt = :central), [Route(:central, sys.A)])[1].kind === :bolus   # bolus-only route is fine
end

@testset "resolve: F, lag, units, event order" begin
    sys = onecmt()
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    # (f) fixed-rate infusion, F = 0.5: same rate, half the duration, half the amount
    ev = resolve(Regimen(100.0; cmt = :central, time = 1.0, rate = 10.0), routes, DoseControl(bioav = Dict(:central => 0.5)))
    @test [(e.time, e.kind, e.rate) for e in ev] == [(1.0, :infusion_start, 10.0), (6.0, :infusion_stop, 10.0)]
    # fixed-duration infusion, F = 0.5: same duration, half the rate
    ev = resolve(Regimen(100.0; cmt = :central, time = 1.0, duration = 10.0), routes, DoseControl(bioav = Dict(:central => 0.5)))
    @test [(e.time, e.kind, e.rate) for e in ev] == [(1.0, :infusion_start, 5.0), (11.0, :infusion_stop, 5.0)]
    # lag shifts the start and the stop; F does not shift anything
    ev = resolve(Regimen(100.0; cmt = :central, time = 1.0, rate = 10.0), routes, DoseControl(lag = Dict(:central => 0.25)))
    @test [e.time for e in ev] == [1.25, 11.25]
    # bolus: F scales the amount, lag the time
    ev = resolve(Regimen(100.0; cmt = :central, time = 3.0), routes, DoseControl(bioav = Dict(:central => 0.8), lag = Dict(:central => 0.5)))
    @test length(ev) == 1 && ev[1].kind === :bolus && ev[1].time == 3.5 && ev[1].amount ≈ 80.0 && ev[1].rate == 0.0
    # F and lag apply per compartment only
    two = [Route(:depot, sys.A), Route(:central, sys.A; rate = sys.Rin)]
    ev = resolve(Regimen(Regimen(10.0; cmt = :depot), Regimen(10.0; cmt = :central)), two, DoseControl(bioav = Dict(:depot => 0.5)))
    @test [(e.cmt, e.amount) for e in ev] == [(:depot, 5.0), (:central, 10.0)]
    # F = 0 on an infusion removes it
    @test isempty(resolve(Regimen(100.0; cmt = :central, rate = 10.0), routes, DoseControl(bioav = Dict(:central => 0.0))))
    # coincident order: stop, start, bolus (input order must not matter)
    r = Regimen(Regimen(5.0; cmt = :central, time = 2.0),
        Regimen(20.0; cmt = :central, time = 2.0, rate = 10.0),      # starts at 2
        Regimen(20.0; cmt = :central, time = 0.0, rate = 10.0))      # stops at 2
    ev = resolve(r, routes)
    @test [e.kind for e in ev if e.time == 2.0] == [:infusion_stop, :infusion_start, :bolus]
    @test issorted(ev; by = e -> e.time)
end

@testset "(b) repeated IV bolus vs closed form" begin
    sys = onecmt()
    routes = [Route(:central, sys.A)]
    reg = Regimen(50.0; cmt = :central, ii = 6.0, addl = 3)
    ts = collect(0.37:0.5:30.0)
    sol, _ = run_regimen(sys, routes, reg, (0.0, 30.0))
    exact(x) = sum(50.0 * exp(-0.3 * (x - td)) for td in 0.0:6.0:18.0 if td <= x; init = 0.0)
    @test all(close(sol(x; idxs = sys.A), exact(x)) for x in ts)
end

@testset "(c) zero-order infusion vs closed form" begin
    sys = onecmt()
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    R, k, V = 8.0, 0.3, 5.0
    sol, prob = run_regimen(sys, routes, Regimen(40.0; cmt = :central, time = 0.0, rate = R), (0.0, 20.0))   # runs 0..5
    xs = collect(0.13:0.41:19.9)
    @test all(close(sol(x; idxs = sys.A) / V, R / (k * V) * (1 - exp(-k * x))) for x in xs if x < 5.0)
    Cend = R / (k * V) * (1 - exp(-k * 5.0))
    @test all(close(sol(x; idxs = sys.A) / V, Cend * exp(-k * (x - 5.0))) for x in xs if x > 5.0)
    @test getp(sys, sys.Rin)(sol.prob) == 0.0           # solved problem is left as it was (sol.prob shares p with the integrator)
    # same result on a second solve of the same problem
    @test solve(prob, Tsit5(); TOL...)(10.0; idxs = sys.A) == sol(10.0; idxs = sys.A)
    # stiff solver gives the same answer
    sol2, _ = run_regimen(sys, routes, Regimen(40.0; cmt = :central, rate = R), (0.0, 20.0); alg = Rodas5P(), tol = (abstol = 1e-11, reltol = 1e-11))
    @test all(close(sol2(x; idxs = sys.A), infusion_exact(x, 0.0, 5.0, R, k); rtol = 1e-6) for x in xs)
end

@testset "(d) overlapping infusions: superposition" begin
    sys = onecmt()
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    k = 0.3
    reg = Regimen(Regimen(30.0; cmt = :central, time = 0.0, rate = 5.0),      # rate 5 on [0, 6]
        Regimen(24.0; cmt = :central, time = 3.0, rate = 4.0),                # rate 4 on [3, 9]
        Regimen(10.0; cmt = :central, time = 4.0, duration = 10.0))           # rate 1 on [4, 14]
    sol, prob = run_regimen(sys, routes, reg, (0.0, 30.0))
    exact(x) = infusion_exact(x, 0.0, 6.0, 5.0, k) + infusion_exact(x, 3.0, 9.0, 4.0, k) + infusion_exact(x, 4.0, 14.0, 1.0, k)
    xs = collect(0.21:0.37:29.9)
    @test all(close(sol(x; idxs = sys.A), exact(x)) for x in xs)
    @test getp(sys, sys.Rin)(sol.prob) == 0.0
    # back-to-back infusions: a stop and a start at the same time give a continuous rate
    back = Regimen(Regimen(6.0; cmt = :central, time = 0.0, rate = 3.0), Regimen(6.0; cmt = :central, time = 2.0, rate = 3.0))
    sol, _ = run_regimen(sys, routes, back, (0.0, 12.0))
    @test all(close(sol(x; idxs = sys.A), infusion_exact(x, 0.0, 4.0, 3.0, k)) for x in 0.3:0.7:11.9)
end

@testset "rate parameter returns to its pre-infusion value" begin
    # k = 0 and a baseline input of 0.5: A(t) = 0.5 t + infused amount
    sys = onecmt(kval = 0.0, rin = 0.5)
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    reg = Regimen(Regimen(10.0; cmt = :central, time = 1.0, rate = 2.0), Regimen(6.0; cmt = :central, time = 3.0, rate = 3.0))   # 2 on [1,6], 3 on [3,5]
    sol, prob = run_regimen(sys, routes, reg, (0.0, 10.0))
    @test close(sol(10.0; idxs = sys.A), 0.5 * 10 + 10.0 + 6.0)
    @test close(sol(4.0; idxs = sys.A), 0.5 * 4 + 2.0 * 3 + 3.0 * 1)
    @test getp(sys, sys.Rin)(sol.prob) == 0.5
    # stopping inside an infusion leaves the problem untouched too
    sol, prob = run_regimen(sys, routes, reg, (0.0, 4.0))
    @test getp(sys, sys.Rin)(sol.prob) == 0.5
    @test close(sol(4.0; idxs = sys.A), 0.5 * 4 + 2.0 * 3 + 3.0 * 1)
end

@testset "(e) oral depot with bioavailability and lag" begin
    ka, k, F, lag, D0 = 1.1, 0.3, 0.7, 0.8, 100.0
    sys = oral(kaval = ka, kval = k)
    routes = [Route(:depot, sys.Ad), Route(:central, sys.Ac; rate = sys.Rin)]
    ctl = DoseControl(bioav = Dict(:depot => F), lag = Dict(:depot => lag))
    sol, _ = run_regimen(sys, routes, Regimen(D0; cmt = :depot, time = 1.0), (0.0, 30.0); control = ctl)
    tl = 1.0 + lag
    bateman(x) = x <= tl ? 0.0 : F * D0 * ka / (ka - k) * (exp(-k * (x - tl)) - exp(-ka * (x - tl)))
    depot(x) = x <= tl ? 0.0 : F * D0 * exp(-ka * (x - tl))
    xs = collect(0.3:0.45:29.8)
    @test all(close(sol(x; idxs = sys.Ac), bateman(x)) for x in xs)
    @test all(close(sol(x; idxs = sys.Ad), depot(x)) for x in xs)
    # repeated oral doses: superposition of Bateman functions
    reg = Regimen(D0; cmt = :depot, ii = 12.0, addl = 2)
    sol, _ = run_regimen(sys, routes, reg, (0.0, 40.0); control = DoseControl(bioav = Dict(:depot => F)))
    multi(x) = sum(x <= td ? 0.0 : F * D0 * ka / (ka - k) * (exp(-k * (x - td)) - exp(-ka * (x - td))) for td in (0.0, 12.0, 24.0))
    @test all(close(sol(x; idxs = sys.Ac), multi(x)) for x in 0.3:0.77:39.9)
end

@testset "(f) F on infusions, checked by amount delivered" begin
    # k = 0: A(t_end) is the total amount delivered
    sys = onecmt(kval = 0.0)
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    F = 0.6
    ctl = DoseControl(bioav = Dict(:central => F))
    # fixed rate 10, amt 100: end time 1 + F*10 = 7, amount 60
    sol, _ = run_regimen(sys, routes, Regimen(100.0; cmt = :central, time = 1.0, rate = 10.0), (0.0, 20.0); control = ctl)
    @test close(sol(20.0; idxs = sys.A), 60.0)
    @test close(sol(6.5; idxs = sys.A), 10.0 * 5.5)       # still infusing at rate 10
    @test close(sol(7.5; idxs = sys.A), 60.0)             # already stopped at 7
    # fixed duration 10, amt 100: ends at 11, rate 6, amount 60
    sol, _ = run_regimen(sys, routes, Regimen(100.0; cmt = :central, time = 1.0, duration = 10.0), (0.0, 20.0); control = ctl)
    @test close(sol(20.0; idxs = sys.A), 60.0)
    @test close(sol(6.0; idxs = sys.A), 6.0 * 5.0)
    @test close(sol(10.5; idxs = sys.A), 6.0 * 9.5) && close(sol(11.5; idxs = sys.A), 60.0)
end

@testset "(g) coincident events all apply" begin
    sys = onecmt(kval = 0.0)
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    # old bolus_callback bug: the second dose at the same time was dropped
    cb = bolus_callback(sys, sys.A, [0.0, 5.0, 5.0, 5.0], [1.0, 2.0, 3.0, 4.0])
    sol = solve(ODEProblem(sys, [], (0.0, 10.0); callback = cb), Tsit5(); TOL...)
    @test close(sol(10.0; idxs = sys.A), 1.0 + 2.0 + 3.0 + 4.0)
    @test close(sol(4.0; idxs = sys.A), 1.0)
    # same through Regimen: two regimens dosing the same compartment at the same time
    reg = Regimen(Regimen(3.0; cmt = :central, time = 2.0), Regimen(4.0; cmt = :central, time = 2.0), Regimen(5.0; cmt = :central, time = 2.0))
    sol, _ = run_regimen(sys, routes, reg, (0.0, 10.0))
    @test close(sol(10.0; idxs = sys.A), 12.0)
    # bolus plus infusion start at the same time
    reg = Regimen(Regimen(7.0; cmt = :central, time = 1.0), Regimen(10.0; cmt = :central, time = 1.0, rate = 2.0))
    sol, _ = run_regimen(sys, routes, reg, (0.0, 10.0))
    @test close(sol(10.0; idxs = sys.A), 17.0)
    # events a rounding error apart stay separate events (only exactly equal times group)
    ev = resolve(Regimen(Regimen(0.2; cmt = :central, time = 0.1, rate = 1.0), Regimen(1.0; cmt = :central, time = 0.3)), routes)
    @test ev[2].time != ev[3].time                      # 0.1 + 0.2 is not 0.3
    cb, tstops = dosing_callback(sys, ev, routes)
    @test length(tstops) == 3
    sol = solve(ODEProblem(sys, [], (0.0, 2.0); callback = cb, tstops = tstops), Tsit5(); TOL...)
    @test close(sol(2.0; idxs = sys.A), 0.2 + 1.0)
end

@testset "(h) initial time and tspan edges" begin
    sys = onecmt()
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    k = 0.3
    # dose at t0 is applied
    sol, _ = run_regimen(sys, routes, Regimen(10.0; cmt = :central, time = 0.0), (0.0, 10.0))
    @test close(sol(2.0; idxs = sys.A), 10.0 * exp(-k * 2.0))
    # infusion at t0 is applied
    sol, _ = run_regimen(sys, routes, Regimen(10.0; cmt = :central, time = 0.0, rate = 5.0), (0.0, 10.0))
    @test close(sol(1.0; idxs = sys.A), infusion_exact(1.0, 0.0, 2.0, 5.0, k))
    # doses after tspan are ignored (bolus and infusion), and nothing breaks
    reg = Regimen(Regimen(10.0; cmt = :central, time = 12.0), Regimen(10.0; cmt = :central, time = 15.0, rate = 1.0),
        Regimen(10.0; cmt = :central, time = 1.0))
    sol, prob = run_regimen(sys, routes, reg, (0.0, 10.0))
    @test close(sol(10.0; idxs = sys.A), 10.0 * exp(-k * 9.0))
    @test sol.retcode == ReturnCode.Success && getp(sys, sys.Rin)(sol.prob) == 0.0
    # a lagged dose pushed past the end is ignored too
    sol, _ = run_regimen(sys, routes, Regimen(10.0; cmt = :central, time = 8.0), (0.0, 10.0); control = DoseControl(lag = Dict(:central => 5.0)))
    @test sol(10.0; idxs = sys.A) == 0.0
    # tspan starting later: a bolus before it is ignored, an infusion still running is on
    reg = Regimen(Regimen(10.0; cmt = :central, time = 0.0), Regimen(40.0; cmt = :central, time = 1.0, rate = 4.0))   # rate 4 on [1, 11]
    sol, prob = run_regimen(sys, routes, reg, (2.0, 20.0))
    @test close(sol(5.0; idxs = sys.A), 4.0 / k * (1 - exp(-k * 3.0)))
    @test close(sol(15.0; idxs = sys.A), 4.0 / k * (1 - exp(-k * 9.0)) * exp(-k * 4.0))
    @test getp(sys, sys.Rin)(sol.prob) == 0.0
    # empty regimen
    sol, _ = run_regimen(sys, routes, Regimen(), (0.0, 5.0))
    @test sol(5.0; idxs = sys.A) == 0.0
end

@testset "(i) route unit scale" begin
    sys = onecmt()
    # dose in mg, model in ug: scale 1000
    routes = [Route(:central, sys.A; rate = sys.Rin, scale = 1000.0)]
    ev = resolve(Regimen(2.0; cmt = :central), routes)
    @test ev[1].amount == 2000.0
    sol, _ = run_regimen(sys, routes, Regimen(2.0; cmt = :central), (0.0, 10.0))
    @test close(sol(3.0; idxs = sys.A), 2000.0 * exp(-0.9))
    # infusion: 10 mg at 5 mg/h is 5000 ug/h for 2 h
    ev = resolve(Regimen(10.0; cmt = :central, rate = 5.0), routes)
    @test [(e.time, e.rate) for e in ev] == [(0.0, 5000.0), (2.0, 5000.0)]
    sol, _ = run_regimen(sys, routes, Regimen(10.0; cmt = :central, rate = 5.0), (0.0, 10.0))
    @test close(sol(6.0; idxs = sys.A), infusion_exact(6.0, 0.0, 2.0, 5000.0, 0.3))
end

@testset "exact restoration and rate arithmetic" begin
    # non-dyadic baseline: 0.2 + 0.1 - 0.1 is not 0.2 in floating point; the solve ends mid-infusion
    sys = onecmt(kval = 0.0, rin = 0.2)
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    sol, _ = run_regimen(sys, routes, Regimen(1.0; cmt = :central, time = 1.0, rate = 0.1), (0.0, 4.0))
    @test getp(sys, sys.Rin)(sol.prob) === 0.2
    sol, _ = run_regimen(sys, routes, Regimen(1.0; cmt = :central, time = 1.0, rate = 0.1), (0.0, 40.0))   # infusion completes
    @test getp(sys, sys.Rin)(sol.prob) === 0.2
    @test close(sol(40.0; idxs = sys.A), 0.2 * 40 + 1.0)

    # mixed scales: 1e16 on [1, 2] and 1 on [0, 10]; after the big one stops the small one is exactly 1
    sys = onecmt(kval = 0.0)
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    seen = Dict{Float64, Float64}()
    probe = DiscreteCallback((u, t_, i) -> t_ in (0.5, 1.5, 3.0), i -> (seen[i.t] = getp(sys, sys.Rin)(i)))
    reg = Regimen(Regimen(1e16; cmt = :central, time = 1.0, rate = 1e16), Regimen(10.0; cmt = :central, time = 0.0, rate = 1.0))
    ev = resolve(reg, routes)
    cb, ts = dosing_callback(sys, ev, routes)
    sol = solve(ODEProblem(sys, [], (0.0, 5.0); callback = CallbackSet(cb, probe), tstops = vcat(ts, 0.5, 1.5, 3.0)), Tsit5(); abstol = 1e-6, reltol = 1e-6)
    @test seen[0.5] == 1.0 && seen[1.5] == 1e16 + 1.0 && seen[3.0] == 1.0
    @test getp(sys, sys.Rin)(sol.prob) === 0.0

    # infusion that stops exactly at t0 is off; one that starts exactly at t0 is on
    sys = onecmt(kval = 0.0)
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    reg = Regimen(Regimen(10.0; cmt = :central, time = 0.0, rate = 5.0), Regimen(6.0; cmt = :central, time = 2.0, rate = 3.0))   # 5 on [0,2], 3 on [2,4]
    sol, _ = run_regimen(sys, routes, reg, (2.0, 10.0))
    @test close(sol(10.0; idxs = sys.A), 6.0) && sol(1.0 + 2.0; idxs = sys.A) > 0
    @test getp(sys, sys.Rin)(sol.prob) === 0.0

    # an infusion much shorter than the roundoff of the time grid still delivers its amount
    reg = Regimen(1.0; cmt = :central, time = 3.0, rate = 1e10)       # lasts 1e-10
    ev = resolve(reg, routes)
    @test ev[2].time > ev[1].time
    sol, _ = run_regimen(sys, routes, reg, (0.0, 6.0))
    @test close(sol(6.0; idxs = sys.A), 1.0; rtol = 1e-6)
    # one so short that start == stop in floating point is rejected, not dropped
    @test_throws ArgumentError resolve(Regimen(1.0; cmt = :central, time = 1e10, rate = 1e10), routes)
    @test_throws ArgumentError dosing_callback(sys, [DoseEvent(1.0, :infusion_stop, :central, 0.0, 2.0)], routes)
    @test_throws ArgumentError dosing_callback(sys, [DoseEvent(1.0, :infusion_start, :central, 0.0, 2.0), DoseEvent(1.0, :infusion_stop, :central, 0.0, 2.0)], routes)
end

@testset "boundaries of tspan" begin
    sys = onecmt(kval = 0.0)
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    amt(t_, a) = Regimen(a; cmt = :central, time = t_)
    # events within 1e-10 of t0 and tend: before t0 and after tend are out, t0 and just inside are in
    reg = Regimen(amt(-5e-10, 100.0), amt(0.0, 1.0), amt(5e-10, 2.0), amt(10.0 - 5e-10, 4.0), amt(10.0 + 5e-10, 1000.0))
    sol, _ = run_regimen(sys, routes, reg, (0.0, 10.0))
    @test close(sol(5.0; idxs = sys.A), 3.0)
    @test sol.u[end][1] == 7.0                       # the dose just inside the end counted
    # a dose exactly at tspan[2] fires (the value after the dose is the last saved point)
    sol, _ = run_regimen(sys, routes, Regimen(amt(2.0, 1.0), amt(10.0, 5.0)), (0.0, 10.0))
    @test sol.u[end][1] == 6.0 && sol.t[end] == 10.0
    # dose exactly at t0 plus one a rounding error before it: only the one at t0 counts
    sol, _ = run_regimen(sys, routes, Regimen(amt(-1e-10, 100.0), amt(0.0, 1.0)), (0.0, 5.0))
    @test close(sol(5.0; idxs = sys.A), 1.0)
    # infusion starting 1e-10 after tspan[2] never turns on; the parameter is untouched
    sol, _ = run_regimen(sys, routes, Regimen(10.0; cmt = :central, time = 10.0 + 1e-10, rate = 1.0), (0.0, 10.0))
    @test sol(10.0; idxs = sys.A) == 0.0 && getp(sys, sys.Rin)(sol.prob) === 0.0
    # infusion stopping just before t0 is off, just after is on at t0
    sol, _ = run_regimen(sys, routes, Regimen(2.0; cmt = :central, time = 0.0, rate = 1.0), (2.0 + 1e-10, 5.0))   # stops at 2
    @test sol(5.0; idxs = sys.A) == 0.0
    sol, _ = run_regimen(sys, routes, Regimen(2.0; cmt = :central, time = 0.0, rate = 1.0), (2.0 - 1e-10, 5.0))
    @test close(sol(5.0; idxs = sys.A), 1e-10; atol = 1e-12)
end

@testset "one callback, several spans and threads" begin
    sys = onecmt()
    routes = [Route(:central, sys.A; rate = sys.Rin)]
    k = 0.3
    reg = Regimen(Regimen(30.0; cmt = :central, time = 1.0, rate = 5.0), Regimen(10.0; cmt = :central, time = 3.0))   # rate 5 on [1,7]; bolus at 3
    ev = resolve(reg, routes)
    cb, ts = dosing_callback(sys, ev, routes)
    exact(x) = infusion_exact(x, 1.0, 7.0, 5.0, k) + (x >= 3.0 ? 10.0 * exp(-k * (x - 3.0)) : 0.0)
    # same callback object, different tspans (before, inside, and after the infusion)
    for tspan in ((0.0, 12.0), (2.0, 12.0), (5.0, 12.0), (7.0, 12.0), (0.0, 4.0), (8.0, 12.0))
        prob = ODEProblem(sys, [sys.A => 0.0], tspan; callback = cb, tstops = ts)
        sol = solve(prob, Tsit5(); TOL...)
        t1 = tspan[1]
        # what the model holds at the start of the span depends on the doses before it, which are not
        # replayed: compare increments of the part that starts inside the span
        if t1 == 0.0
            @test all(close(sol(x; idxs = sys.A), exact(x)) for x in (tspan[2] == 4.0 ? (0.5, 1.7, 3.3, 3.9) : (0.5, 2.2, 3.5, 6.1, 9.4, 11.7)))
        end
        @test getp(sys, sys.Rin)(sol.prob) === 0.0
    end
    # infusion already running at t0 = 2 starts the model from 0 with the input on
    sol = solve(ODEProblem(sys, [sys.A => 0.0], (2.0, 12.0); callback = cb, tstops = ts), Tsit5(); TOL...)
    inf(x) = infusion_exact(x, 2.0, 7.0, 5.0, k)
    @test close(sol(2.5; idxs = sys.A), inf(2.5)) && close(sol(6.0; idxs = sys.A), inf(6.0) + 10.0 * exp(-k * 3.0))

    # threads: ONE callback shared by concurrent solves with different k and baseline input
    ks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    bases = [0.0, 0.5, 0.25, 1.0, 0.0, 0.7, 0.2, 0.9]
    base = ODEProblem(sys, [], (0.0, 12.0); callback = cb, tstops = ts)
    ens = EnsembleProblem(base; prob_func = (prob, ctx) -> remake(prob; p = [sys.k => ks[ctx.sim_id], sys.Rin => bases[ctx.sim_id]]),
        output_func = (sol, i) -> ((sol(9.0; idxs = sys.A), getp(sys, sys.Rin)(sol.prob)), false), safetycopy = false)
    out = solve(ens, Tsit5(), EnsembleThreads(); trajectories = length(ks), abstol = 1e-12, reltol = 1e-12)
    for i in eachindex(ks)
        kk, b = ks[i], bases[i]
        want = b / kk * (1 - exp(-kk * 9.0)) + infusion_exact(9.0, 1.0, 7.0, 5.0, kk) + 10.0 * exp(-kk * 6.0)
        @test close(out.u[i][1], want; rtol = 1e-7)
        @test out.u[i][2] === b
    end
end

@testset "bolus_callback keeps its behaviour" begin
    sys = onecmt()
    times = repeated_doses(0.0, 6.0, 3)
    cb = bolus_callback(sys, sys.A, times, 10.0)
    sol = solve(ODEProblem(sys, [], (0.0, 24.0); callback = cb), Tsit5(); TOL...)
    exact(x) = sum(10.0 * exp(-0.3 * (x - td)) for td in times if td <= x; init = 0.0)
    @test all(close(sol(x; idxs = sys.A), exact(x)) for x in 0.5:1.07:23.5)
    # zero-amount doses only add stop times
    cb0 = bolus_callback(sys, sys.A, [2.5, 7.5], 0.0)
    sol0 = solve(ODEProblem(sys, [sys.A => 1.0], (0.0, 10.0); callback = cb0), Tsit5(); TOL...)
    @test 2.5 in sol0.t && 7.5 in sol0.t && close(sol0(10.0; idxs = sys.A), exp(-3.0))
    @test_throws DimensionMismatch bolus_callback(sys, sys.A, [1.0, 2.0], [1.0])
end


@testset "routes may not share a rate parameter" begin
    @parameters kel = 0.1 Rin = 0.0
    @variables A(t) = 0.0 B(t) = 0.0
    @mtkcompile shsys = System([D(A) ~ -kel * A + Rin, D(B) ~ -kel * B + Rin], t)
    rts = [Route(:a, shsys.A; rate = shsys.Rin), Route(:b, shsys.B; rate = shsys.Rin)]
    ev = resolve(Regimen(1.0; cmt = :a, rate = 1.0), rts)
    @test_throws ArgumentError dosing_callback(shsys, ev, rts)
end
end # module RegimenTests
