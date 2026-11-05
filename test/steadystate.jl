# Steady-state dosing (SS=1) and reset events (EVID=3 / 4). Own module, so helper names cannot
# clash with other test files. Closed forms are the textbook steady-state formulas and, for
# lag / overlapping infusions, a superposition of the exact single-dose responses of the linear
# model (independent of the cycle iteration the package uses).
module SteadyStateTests

using QSPAtlas, Test
using Logging: with_logger, NullLogger
using LinearAlgebra: I, exp
using QSPAtlas: Reset, reset_times, SteadyStateOptions, SteadyStateError, steady_state_cycle, with_steady_state
using ModelingToolkit, OrdinaryDiffEq, SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D

const TOL = (abstol = 1e-11, reltol = 1e-11)

# depot -> central, infusion input Rin into central; k = CL / V
@parameters ka = 1.1 CL = 2.0 V = 10.0 Rin = 0.0
@variables Ad(t) = 0.0 Ac(t) = 0.0
@mtkcompile pk = System([D(Ad) ~ -ka * Ad, D(Ac) ~ ka * Ad - CL / V * Ac + Rin], t)

const ROUTES = [Route(:depot, pk.Ad), Route(:central, pk.Ac; rate = pk.Rin)]

function pkmodel(; Ω = [0.09;;], ηnames = [:CL], init = (θ, η, c) -> Pair[], tspan = (0.0, 100.0), solve_kwargs = TOL)
    return PopulationModel(pk;
        routes = ROUTES,
        θ = (CL = 2.0, V = 10.0, ka = 1.1),
        ηnames, Ω,
        params = (θ, η, c) -> [pk.CL => θ.CL * exp(get(η, :CL, 0.0)) * (get(c, :WT, 70.0) / 70)^0.75, pk.V => θ.V, pk.ka => θ.ka],
        init,
        observed = (A = pk.Ac, Dep = pk.Ad),
        error = (A = AdditiveError(0.1), Dep = AdditiveError(0.1)),
        tspan, solver = Tsit5(), solve_kwargs)
end

# A deterministic model (no random effects) for an arbitrary system.
function simplemodel(sys, routes, pp = Pair[]; observed, solve_kwargs = TOL, tspan = (0.0, 100.0))
    return PopulationModel(sys; routes, θ = (;), ηnames = Symbol[], Ω = zeros(0, 0), params = (θ, η, c) -> pp,
        observed, error = map(_ -> AdditiveError(0.1), observed), tspan, solver = Tsit5(), solve_kwargs)
end

const K = 0.2          # CL / V at eta = 0
const KA = 1.1

# --- exact responses of ONE dose given at time 0, at time x (zero before the dose lands) ---
bolus_c(x; amt, lag = 0.0) = x < lag ? 0.0 : amt * exp(-K * (x - lag))
oral_c(x; amt, lag = 0.0, k = K) = x < lag ? 0.0 : amt * KA / (KA - k) * (exp(-k * (x - lag)) - exp(-KA * (x - lag)))
oral_d(x; amt, lag = 0.0) = x < lag ? 0.0 : amt * exp(-KA * (x - lag))
function infusion_c(x; amt, dur, lag = 0.0, k = K)
    R = amt / dur
    x <= lag && return 0.0
    x <= lag + dur && return R / k * (1 - exp(-k * (x - lag)))
    return R / k * (1 - exp(-k * dur)) * exp(-k * (x - lag - dur))
end

# Superposition: doses at m * tau for m = -M ... addl, observed s after the steady-state dose.
ss_sum(resp, tau, s; addl = 0, M = 300) = sum(resp(s - m * tau) for m in -M:addl)

rows(res, endpoint) = [r.ipred for r in res.obs if r.endpoint === endpoint]
const OBS_S = [0.0, 0.7, 2.0, 3.0, 5.0, 5.999, 6.0, 8.5, 12.0, 17.0]      # s after the steady-state dose

# One steady-state case: cmt/regimen keywords + control; checks central (and depot) against
# the superposition exact solution. T is the steady-state time (> 0 exercises the jump).
function check_ss(; amt, cmt, tau, T = 0.0, addl = 0, lag = 0.0, F = 1.0, rate = 0.0, duration = 0.0,
        resp_c, resp_d = nothing, rtol = 1e-6, obs_s = OBS_S, label = "", model = pkmodel(), offset = 0.0, atol = 1e-8)
    reg = Regimen(amt; cmt, time = T, ii = tau, addl, ss = true, rate, duration)
    control = DoseControl(bioav = Dict(cmt => F), lag = Dict(cmt => lag))
    s = Subject("s"; regimen = reg, obstimes = T .+ obs_s, control)
    res = zero_re_prediction(model, s)
    got = rows(res, :A)
    want = [offset + ss_sum(resp_c, tau, x; addl) for x in obs_s]
    ok = all(isapprox.(got, want; rtol, atol))
    ok || @info "mismatch $label" got want
    @test ok
    if resp_d !== nothing
        gd = rows(res, :Dep)
        wd = [ss_sum(resp_d, tau, x; addl) for x in obs_s]
        @test all(isapprox.(gd, wd; rtol, atol))
    end
    return res
end

@testset "regimen record: ss and reset" begin
    r = Regimen(100.0; cmt = :central, ii = 12.0, addl = 2, ss = true)
    @test [d.ss for d in doses(r)] == [true, false, false]          # only the first dose is steady state
    @test all(d -> d.ii == 12.0, doses(r)[1:1]) && doses(r)[2].ii == 0.0
    @test_throws ArgumentError Regimen(100.0; cmt = :central, ss = true)                       # ii needed
    @test_throws ArgumentError Regimen(0.0; cmt = :central, ii = 12.0, ss = true)             # amt > 0 needed
    @test_throws ArgumentError Regimen(100.0; cmt = :central, ii = 12.0, ss = true, reset = true)
    @test_throws ArgumentError Dose(0.0, 1.0; cmt = :a, ss = true)
    # two steady-state records with different ii or time cannot be combined; the same ii and time can
    a = Regimen(100.0; cmt = :central, ii = 12.0, ss = true)
    @test_throws ArgumentError vcat(a, Regimen(5.0; cmt = :depot, ii = 24.0, ss = true))
    @test_throws ArgumentError vcat(a, Regimen(5.0; cmt = :depot, ii = 12.0, ss = true, time = 1.0))
    @test length(doses(vcat(a, Regimen(5.0; cmt = :depot, ii = 12.0, ss = true)))) == 2
    # resets: pure, with a dose, combined, de-duplicated, sorted
    @test reset_times(Regimen(Reset(5.0))) == [5.0]
    @test isempty(doses(Regimen(Reset(5.0))))
    rr = Regimen(Regimen(10.0; cmt = :depot, ii = 4.0, addl = 3), Reset(9.0), Regimen(1.0; cmt = :depot, time = 2.0, reset = true))
    @test reset_times(rr) == [2.0, 9.0]
    @test reset_times(vcat(Regimen(Reset(3.0)), Regimen(Reset(3.0)))) == [3.0]
    @test_throws ArgumentError Reset(NaN)
    @test reset_times(Regimen(1.0; cmt = :depot)) == Float64[]
end

@testset "resolve: reset cancels pending doses and running infusions" begin
    ev(reg, ctl = DoseControl()) = [(e.time, e.kind, e.cmt, e.amount, e.rate) for e in resolve(reg, ROUTES, ctl)]
    # infusion rate 10 from 0 (duration 10), reset at 4: stopped at 4, start kept
    reg = Regimen(Regimen(100.0; cmt = :central, duration = 10.0), Reset(4.0))
    @test ev(reg) == [(0.0, :infusion_start, :central, 0.0, 10.0), (4.0, :reset, :all, 0.0, 0.0), (4.0, :infusion_stop, :central, 0.0, 10.0)]
    # an infusion that starts at or after the reset (through lag) is cancelled, one given at the reset time is kept
    ctl = DoseControl(lag = Dict(:central => 3.0))
    e1 = ev(Regimen(Regimen(100.0; cmt = :central, time = 2.0, duration = 10.0), Reset(5.0)), ctl)     # starts 5.0 = reset: cancelled
    @test e1 == [(5.0, :reset, :all, 0.0, 0.0)]
    e2 = ev(Regimen(Regimen(100.0; cmt = :central, time = 5.0, duration = 10.0), Reset(5.0)), ctl)     # given at the reset: kept, lag after it
    @test e2 == [(5.0, :reset, :all, 0.0, 0.0), (8.0, :infusion_start, :central, 0.0, 10.0), (18.0, :infusion_stop, :central, 0.0, 10.0)]
    # boluses: pending at the reset (lands at or after) cancelled; landed before it, kept
    ctl = DoseControl(lag = Dict(:depot => 3.0))
    d = Regimen(5.0; cmt = :depot, ii = 1.0, addl = 3)             # given at 0,1,2,3 land at 3,4,5,6
    e3 = ev(Regimen(d, Reset(4.0)), ctl)
    @test [(x[1], x[2]) for x in e3] == [(3.0, :bolus), (4.0, :reset)]
    # reset + dose (EVID=4): the dose at the reset time follows the reset
    e4 = ev(Regimen(Regimen(7.0; cmt = :depot, time = 4.0, reset = true)))
    @test e4 == [(4.0, :reset, :all, 0.0, 0.0), (4.0, :bolus, :depot, 7.0, 0.0)]
    # a later reset also cancels what an earlier reset left (each dose is cut at the first reset after it was given)
    e5 = ev(Regimen(Regimen(100.0; cmt = :central, duration = 10.0), Reset(2.0), Reset(6.0)))
    @test [(x[1], x[2]) for x in e5] == [(0.0, :infusion_start), (2.0, :reset), (2.0, :infusion_stop), (6.0, :reset)]
end

@testset "resolve: steady-state carry-over events" begin
    ev(reg, ctl = DoseControl()) = [(e.time, e.kind, e.cmt, e.amount, e.rate) for e in resolve(reg, ROUTES, ctl)]
    # bolus with lag 8 > ii = 6: the dose of one interval earlier lands at 2 (still pending at T = 0)
    reg = Regimen(10.0; cmt = :depot, ii = 6.0, ss = true)
    ctl = DoseControl(lag = Dict(:depot => 8.0))
    @test ev(reg, ctl) == [(0.0, :steady_state, :all, 0.0, 0.0), (2.0, :bolus, :depot, 10.0, 0.0), (8.0, :bolus, :depot, 10.0, 0.0)]
    T, ii, cyc = steady_state_cycle(reg, ROUTES, ctl)
    @test (T, ii) == (0.0, 6.0)
    @test [(e.time, e.kind) for e in cyc] == [(2.0, :bolus)]              # the cycle sees one landing per interval
    # infusion of 10 h > ii = 6 h: the dose before is still running at T and ends at 4
    reg = Regimen(100.0; cmt = :central, ii = 6.0, duration = 10.0, ss = true)
    @test [(x[1], x[2]) for x in ev(reg)] == [(0.0, :steady_state), (0.0, :infusion_start), (0.0, :infusion_start), (4.0, :infusion_stop), (10.0, :infusion_stop)]
    _, _, cyc = steady_state_cycle(reg, ROUTES)
    @test [(e.time, e.kind) for e in cyc] == [(0.0, :infusion_start), (0.0, :infusion_start), (4.0, :infusion_stop), (10.0, :infusion_stop)]
    # infusion from a lagged earlier dose is switched on at T (not at its original start)
    reg = Regimen(100.0; cmt = :central, ii = 6.0, duration = 2.0, ss = true, time = 20.0)
    e = ev(reg, DoseControl(lag = Dict(:central => 5.0)))                # dose at 14 starts 19, stops 21; dose at 20 starts 25
    @test [(x[1], x[2]) for x in e] == [(20.0, :steady_state), (20.0, :infusion_start), (21.0, :infusion_stop), (25.0, :infusion_start), (27.0, :infusion_stop)]
    # no steady state: no cycle
    @test steady_state_cycle(Regimen(1.0; cmt = :depot), ROUTES) === nothing
    # a later reset cancels the carried-over events too
    reg = Regimen(Regimen(10.0; cmt = :depot, ii = 6.0, ss = true), Reset(1.0))
    @test [(x[1], x[2]) for x in ev(reg, DoseControl(lag = Dict(:depot => 8.0)))] == [(0.0, :steady_state), (1.0, :reset)]
    # a reset at the steady-state time is superseded
    reg = Regimen(Regimen(10.0; cmt = :depot, ii = 6.0, ss = true), Reset(0.0))
    @test [(x[1], x[2]) for x in ev(reg)] == [(0.0, :steady_state), (0.0, :bolus)]
    # a steady-state event needs its state
    @test_throws ArgumentError dosing_callback(pk, resolve(Regimen(10.0; cmt = :depot, ii = 6.0, ss = true), ROUTES), ROUTES)
    evs = with_steady_state(resolve(Regimen(10.0; cmt = :depot, ii = 6.0, ss = true), ROUTES), [1.0, 2.0])
    @test dosing_callback(pk, evs, ROUTES) isa Tuple
    @test_throws ArgumentError dosing_callback(pk, with_steady_state(evs, [1.0]), ROUTES)       # wrong length
    both = vcat(evs, QSPAtlas.DoseEvent(0.0, :reset, :all, 0.0, 0.0))
    @test_throws ArgumentError dosing_callback(pk, both, ROUTES)
end

@testset "IV bolus at steady state: peak and trough closed forms" begin
    tau, amt = 6.0, 100.0
    peak = amt / (1 - exp(-K * tau))
    trough = peak * exp(-K * tau)
    s = Subject("a"; regimen = Regimen(amt; cmt = :central, ii = tau, ss = true, addl = 2), obstimes = [0.0, tau - 1e-6, tau, 2tau, 3tau])
    r = zero_re_prediction(pkmodel(), s)
    a = rows(r, :A)
    @test a[1] ≈ peak rtol = 1e-7                       # right after the steady-state dose
    @test a[2] ≈ trough rtol = 1e-6                     # just before the next dose
    @test a[3] ≈ peak rtol = 1e-7
    @test a[4] ≈ peak rtol = 1e-7
    @test a[5] ≈ peak * exp(-K * tau) rtol = 1e-7        # the addl doses end; one more interval of decay
    # the same subject is not affected by its initial condition
    m50 = pkmodel(init = (θ, η, c) -> [pk.Ac => 50.0, pk.Ad => 3.0])
    @test rows(zero_re_prediction(m50, s), :A)[1:4] ≈ a[1:4] rtol = 1e-6
    # steady-state dose later in the run: the jump at T supersedes what was given before
    T = 30.0
    reg = Regimen(Regimen(500.0; cmt = :central, time = 0.0), Regimen(amt; cmt = :central, time = T, ii = tau, ss = true))
    s2 = Subject("b"; regimen = reg, obstimes = [29.0, 29.999, T, T + 3.0, T + tau - 1e-6])
    a2 = rows(zero_re_prediction(pkmodel(), s2), :A)
    @test a2[1] ≈ 500 * exp(-K * 29) rtol = 1e-7        # before T: the first dose only
    @test a2[3] ≈ peak rtol = 1e-7                      # at T the state is the steady state, not 500 e^{-kT} + peak
    @test a2[4] ≈ peak * exp(-K * 3) rtol = 1e-7
    @test a2[5] ≈ trough rtol = 1e-6
end

@testset "SS=1 cancels pending ADDL doses of earlier records (NONMEM)" begin
    ev(reg, ctl = DoseControl()) = [(e.time, e.kind) for e in resolve(reg, ROUTES, ctl)]
    # bolus 100 every 12 h from 0 (ADDL 5: 0, 12, ..., 60); SS=1 dose at 30 every 24 h with one ADDL (54)
    early = Regimen(100.0; cmt = :central, ii = 12.0, addl = 5)
    ss = Regimen(50.0; cmt = :central, time = 30.0, ii = 24.0, addl = 1, ss = true)
    @test ev(vcat(early, ss)) == [(0.0, :bolus), (12.0, :bolus), (24.0, :bolus), (30.0, :steady_state), (30.0, :bolus), (54.0, :bolus)]
    # an ADDL dose at the SS time itself is pending too; a later record (and its ADDL) is kept
    later = Regimen(7.0; cmt = :central, time = 40.0, ii = 12.0, addl = 1)
    @test ev(vcat(Regimen(100.0; cmt = :central, ii = 15.0, addl = 3), ss, later)) ==
          [(0.0, :bolus), (15.0, :bolus), (30.0, :steady_state), (30.0, :bolus), (40.0, :bolus), (52.0, :bolus), (54.0, :bolus)]
    # a dose record of its own at a later time (not an ADDL dose of an earlier record) is kept
    @test ev(vcat(Regimen(100.0; cmt = :central, time = 36.0), ss)) == [(30.0, :steady_state), (30.0, :bolus), (36.0, :bolus), (54.0, :bolus)]
    # ADDL infusions: the one running at 30 is stopped there, the later ones are cancelled
    inf = Regimen(100.0; cmt = :central, ii = 12.0, addl = 3, duration = 10.0)    # 0-10, 12-22, 24-34, 36-46
    @test ev(vcat(inf, ss))[end - 4:end] == [(24.0, :infusion_start), (30.0, :steady_state), (30.0, :infusion_stop), (30.0, :bolus), (54.0, :bolus)]
    # an SS=2 record cancels nothing: all six ADDL boluses plus its own dose
    ss2 = Regimen(50.0; cmt = :central, time = 30.0, ii = 24.0, ss = 2)
    @test count(e -> e[2] === :bolus, ev(vcat(early, ss2))) == 7
    # closed form after the steady-state record: the steady state of 50 every 24 h alone
    peak = 50 / (1 - exp(-K * 24))
    obst = [29.0, 30.0, 37.0, 49.0, 53.9, 54.0, 61.0]
    sim(reg) = rows(zero_re_prediction(pkmodel(), Subject("c"; regimen = reg, obstimes = obst)), :A)
    got = sim(vcat(early, ss))
    @test got[1] ≈ 100 * sum(exp(-K * (29 - 12j)) for j in 0:2) rtol = 1e-8
    @test got[2:5] ≈ peak .* exp.(-K .* (obst[2:5] .- 30)) rtol = 1e-7
    @test got[6:7] ≈ (peak * exp(-K * 24) + 50) .* exp.(-K .* (obst[6:7] .- 54)) rtol = 1e-7
    # the same as the regimen with the cancelled doses left out by hand
    @test got ≈ sim(vcat(Regimen(100.0; cmt = :central, ii = 12.0, addl = 2), ss)) rtol = 1e-10
end

@testset "oral first order absorption at steady state (with F)" begin
    tau, amt, F = 12.0, 200.0, 0.8
    # closed form: A_c(s) = F D ka/(ka-k) (e^{-ks}/(1-e^{-k tau}) - e^{-ka s}/(1-e^{-ka tau}))
    cf(s) = F * amt * KA / (KA - K) * (exp(-K * s) / (1 - exp(-K * tau)) - exp(-KA * s) / (1 - exp(-KA * tau)))
    df(s) = F * amt * exp(-KA * s) / (1 - exp(-KA * tau))
    ss = [0.0, 1.0, 2.5, 6.0, 11.9, 12.0]
    s = Subject("o"; regimen = Regimen(amt; cmt = :depot, ii = tau, ss = true, addl = 1), obstimes = ss, control = DoseControl(bioav = Dict(:depot => F)))
    r = zero_re_prediction(pkmodel(), s)
    @test rows(r, :A)[1:5] ≈ cf.(ss[1:5]) rtol = 1e-6
    @test rows(r, :Dep)[1:5] ≈ df.(ss[1:5]) rtol = 1e-6
    # the same through the generic superposition check, with a steady-state time away from zero
    check_ss(; amt, cmt = :depot, tau, T = 24.0, F, resp_c = x -> oral_c(x; amt = F * amt), resp_d = x -> oral_d(x; amt = F * amt), label = "oral T=24")
end

@testset "zero-order infusion at steady state: shorter and longer than tau" begin
    k, tau, amt = K, 6.0, 100.0
    # duration d = 2 < tau: trough and in-infusion curve in closed form
    d = 2.0
    R = amt / d
    Aend = R / k * (1 - exp(-k * d)) / (1 - exp(-k * tau))     # amount at the end of the infusion
    Atr = Aend * exp(-k * (tau - d))                           # trough, just before the next start
    s = Subject("i"; regimen = Regimen(amt; cmt = :central, ii = tau, ss = true, duration = d, addl = 1), obstimes = [0.0, 1.0, d, 4.0, tau - 1e-6])
    a = rows(zero_re_prediction(pkmodel(), s), :A)
    @test a[1] ≈ Atr rtol = 1e-6
    @test a[2] ≈ R / k * (1 - exp(-k)) + Atr * exp(-k) rtol = 1e-6
    @test a[3] ≈ Aend rtol = 1e-6
    @test a[4] ≈ Aend * exp(-k * 2) rtol = 1e-6
    @test a[5] ≈ Atr rtol = 1e-5
    # duration 2 by rate instead of duration, and checks by superposition
    check_ss(; amt, cmt = :central, tau, rate = 50.0, resp_c = x -> infusion_c(x; amt, dur = 2.0), label = "inf d<tau (rate)")
    # duration 10 > tau = 6: infusions overlap; the infusion that began one interval earlier is running at T
    check_ss(; amt, cmt = :central, tau, duration = 10.0, addl = 3, resp_c = x -> infusion_c(x; amt, dur = 10.0), label = "inf d>tau")
    # duration 14 > 2 tau (three infusions overlap at times)
    check_ss(; amt, cmt = :central, tau, duration = 14.0, addl = 2, resp_c = x -> infusion_c(x; amt, dur = 14.0), label = "inf d>2tau")
    # overlapped infusion, steady-state dose at T > 0
    check_ss(; amt, cmt = :central, tau, T = 18.0, duration = 10.0, addl = 3, resp_c = x -> infusion_c(x; amt, dur = 10.0), label = "inf d>tau T")
    # F scales the DURATION of a fixed-rate infusion: rate 50, amt 100, F = 0.5 -> 1 h
    check_ss(; amt, cmt = :central, tau, rate = 50.0, F = 0.5, resp_c = x -> infusion_c(x; amt = 50.0, dur = 1.0), label = "inf F")
end

@testset "steady state with lag" begin
    tau, amt = 6.0, 100.0
    # infusion (duration 2) delayed by 3 h: lag + duration = 5 < tau, everything inside the interval
    check_ss(; amt, cmt = :central, tau, lag = 3.0, duration = 2.0, resp_c = x -> infusion_c(x; amt, dur = 2.0, lag = 3.0), label = "inf lag<tau")
    # lag 5 + duration 3: the infusion crosses the interval boundary
    check_ss(; amt, cmt = :central, tau, lag = 5.0, duration = 3.0, addl = 1, resp_c = x -> infusion_c(x; amt, dur = 3.0, lag = 5.0), label = "inf crosses boundary")
    # lag 8 > tau: the whole infusion lands in the next interval, so the dose of the interval before is running at T
    check_ss(; amt, cmt = :central, tau, lag = 8.0, duration = 2.0, addl = 2, resp_c = x -> infusion_c(x; amt, dur = 2.0, lag = 8.0), label = "inf lag>tau")
    # lag 13 > 2 tau, steady-state time away from zero
    check_ss(; amt, cmt = :central, tau, T = 12.0, lag = 13.0, duration = 2.5, addl = 2, resp_c = x -> infusion_c(x; amt, dur = 2.5, lag = 13.0), label = "inf lag>2tau")
    # oral bolus: lag 1 (inside), lag 7 > tau (wraps), lag 12.5 > 2 tau
    for lag in (1.0, 7.0, 12.5)
        check_ss(; amt, cmt = :depot, tau, lag, addl = 1, F = 0.7, resp_c = x -> oral_c(x; amt = 0.7 * amt, lag),
            resp_d = x -> oral_d(x; amt = 0.7 * amt, lag), label = "oral lag=$lag")
    end
    # IV bolus with lag 9
    check_ss(; amt, cmt = :central, tau, lag = 9.0, addl = 1, resp_c = x -> bolus_c(x; amt, lag = 9.0), label = "bolus lag")
end

@testset "reset: states back to initial conditions, infusions and pending doses cancelled" begin
    m = pkmodel(init = (θ, η, c) -> [pk.Ac => 50.0, pk.Ad => 0.0])
    sub(reg; obs, ctl = DoseControl()) = zero_re_prediction(m, Subject("r"; regimen = reg, obstimes = obs, control = ctl))
    # bolus at 0, reset at 5: back to the initial amount 50 (not zero), then plain decay
    r = sub(Regimen(Regimen(100.0; cmt = :central), Reset(5.0)); obs = [1.0, 4.999, 5.0, 7.0])
    @test rows(r, :A) ≈ [150 * exp(-K * 1), 150 * exp(-K * 4.999), 50.0, 50 * exp(-K * 2)] rtol = 1e-7
    # running infusion (rate 10 for 10 h) cancelled by the reset at 4: Rin back to its baseline 0
    r = sub(Regimen(Regimen(100.0; cmt = :central, duration = 10.0), Reset(4.0)); obs = [3.0, 4.0, 6.0, 12.0])
    inf(x) = 50 * exp(-K * x) + 10 / K * (1 - exp(-K * x))
    @test rows(r, :A) ≈ [inf(3.0), 50.0, 50 * exp(-K * 2), 50 * exp(-K * 8)] rtol = 1e-7
    # pending lagged doses: lag 3 on the depot. dose at 0 lands at 3 (before the reset at 4, kept);
    # dose at 2 would land at 5 (cancelled by the reset); dose at 4 (given at the reset time) lands at 7 (kept)
    ctl = DoseControl(lag = Dict(:depot => 3.0))
    reg = Regimen(Regimen(10.0; cmt = :depot, time = 0.0), Regimen(20.0; cmt = :depot, time = 2.0), Reset(4.0), Regimen(40.0; cmt = :depot, time = 4.0))
    r = sub(reg; obs = [3.5, 4.0, 5.0, 6.0, 7.0, 9.0], ctl)
    ad(x) = x < 7 ? 0.0 : 40 * exp(-KA * (x - 7))
    @test rows(r, :Dep) ≈ [10 * exp(-KA * 0.5), 0.0, 0.0, 0.0, 40.0, ad(9.0)] rtol = 1e-7 atol = 1e-10
    # amount in the central compartment: initial 50 after the reset, plus the absorbed dose of 40
    ac(x) = 50 * exp(-K * (x - 4)) + (x < 7 ? 0.0 : 40 * KA / (KA - K) * (exp(-K * (x - 7)) - exp(-KA * (x - 7))))
    @test rows(r, :A)[3:6] ≈ ac.([5.0, 6.0, 7.0, 9.0]) rtol = 1e-6
    # reset + dose (EVID=4): central dose at the reset time sees the initial state plus the dose
    reg = Regimen(Regimen(100.0; cmt = :central), Regimen(30.0; cmt = :central, time = 4.0, reset = true))
    r = sub(reg; obs = [3.0, 4.0, 6.0])
    @test rows(r, :A) ≈ [150 * exp(-K * 3), 80.0, 80 * exp(-K * 2)] rtol = 1e-7
    # reset after a steady-state dose goes back to the initial condition, not to the steady state
    reg = Regimen(Regimen(100.0; cmt = :central, ii = 6.0, ss = true), Reset(9.0))
    r = sub(reg; obs = [0.0, 8.999, 9.0, 12.0])
    peak = 100 / (1 - exp(-K * 6))
    @test rows(r, :A) ≈ [peak, peak * exp(-K * 8.999), 50.0, 50 * exp(-K * 3)] rtol = 1e-6
    # reset at the time of an infusion start: the infusion starts after the reset
    reg = Regimen(Regimen(100.0; cmt = :central, duration = 10.0), Regimen(100.0; cmt = :central, duration = 10.0, time = 3.0, reset = true))
    r = sub(reg; obs = [3.0, 5.0])
    @test rows(r, :A) ≈ [50.0, 50 * exp(-K * 2) + 10 / K * (1 - exp(-K * 2))] rtol = 1e-7
    # a reset after the end of the solve is harmless; an observation at the reset time sees the state after it
    r = sub(Regimen(Regimen(100.0; cmt = :central), Reset(500.0)); obs = [2.0])
    @test rows(r, :A)[1] ≈ 150 * exp(-K * 2) rtol = 1e-7
    r = sub(Regimen(Regimen(100.0; cmt = :central), Reset(10.0)); obs = [10.0])
    @test rows(r, :A)[1] ≈ 50.0 rtol = 1e-9
end

@testset "reset with an infusion leaves the rate parameter at its baseline" begin
    # nonzero baseline rate: 0.5 per hour keeps flowing after the infusion is cancelled
    @parameters Rb = 0.5
    @variables X(t) = 0.0
    @mtkcompile sysb = System([D(X) ~ -0.1 * X + Rb], t)
    routes = [Route(:c, sysb.X; rate = sysb.Rb)]
    reg = Regimen(Regimen(20.0; cmt = :c, duration = 10.0), Reset(5.0))      # extra rate 2
    ev = resolve(reg, routes)
    cb, ts = dosing_callback(sysb, ev, routes)
    prob = ODEProblem(sysb, [], (0.0, 15.0); callback = cb, tstops = ts)
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)
    # 0..5: rate 2.5; reset to X(0) = 0 at 5; 5..15: rate 0.5
    x5 = 2.5 / 0.1 * (1 - exp(-0.1 * 5))
    @test sol(5.0 - 1e-9)[1] ≈ x5 rtol = 1e-6 atol = 1e-9
    @test sol.u[end][1] ≈ 0.5 / 0.1 * (1 - exp(-0.1 * 10)) rtol = 1e-7
end

@testset "steady state: failures are explicit" begin
    @variables B(t) = 0.0
    @parameters g = 0.1
    @mtkcompile grow = System([D(B) ~ g * B], t)
    mg = PopulationModel(grow; routes = [Route(:c, grow.B)], θ = (g = 0.1,), ηnames = Symbol[], Ω = zeros(0, 0),
        params = (θ, η, c) -> [grow.g => θ.g], observed = (B = grow.B,), error = (B = AdditiveError(1.0),),
        tspan = (0.0, 48.0), solver = Tsit5())
    sg = Subject("g"; regimen = Regimen(1.0; cmt = :c, ii = 6.0, ss = true), obstimes = [1.0])
    err = try
        simulate(mg, sg; η = :zero, steady_state = SteadyStateOptions(maxcycles = 50)); nothing
    catch e
        e
    end
    @test err isa SteadyStateError
    @test occursin("no steady state reached after", sprint(showerror, err))
    # in a cohort the failure is a trajectory outcome: flagged rows, NaN observations, one warning,
    # and the other subjects are solved normally
    mgc = PopulationModel(grow; routes = [Route(:c, grow.B)], θ = (g = 0.1,), ηnames = Symbol[], Ω = zeros(0, 0),
        params = (θ, η, c) -> [grow.g => get(c, :G, θ.g)], observed = (B = grow.B,), error = (B = AdditiveError(1.0),),
        tspan = (0.0, 48.0), solver = Tsit5(), solve_kwargs = TOL)
    mk(id, G) = Subject(id; regimen = Regimen(1.0; cmt = :c, ii = 6.0, ss = true), covariates = (G = G,), obstimes = [1.0, 3.0])
    cohort = [mk("ok1", -0.3), mk("bad", 0.1), mk("ok2", -0.2)]
    opts = SteadyStateOptions(maxcycles = 50)
    res = @test_logs (:warn, r"1 of 3 subject simulations") simulate(mgc, cohort; steady_state = opts, residual = false)
    @test [r.retcode for r in res.subjects] == [:Success, :SteadyStateFailed, :Success]
    @test occursin("no steady state reached after 50 cycles", res.subjects[2].message)
    @test !haskey(res.subjects[1], :message)
    bad = [o for o in res.obs if o.id == "bad"]
    @test length(bad) == 2 && all(o -> isnan(o.ipred) && o.retcode === :SteadyStateFailed, bad)
    for (id, g) in (("ok1", -0.3), ("ok2", -0.2))
        got = [o.ipred for o in res.obs if o.id == id]
        @test got ≈ [exp(g * x) / (1 - exp(g * 6.0)) for x in (1.0, 3.0)] rtol = 1e-6
    end
    # replicates, residual draws and serial/threaded execution keep working with a failed trajectory
    r2 = @test_logs (:warn, r"2 of 6 subject simulations") simulate(mgc, cohort; steady_state = opts, replicates = 2, seed = 5)
    r2s = @test_logs (:warn, r"2 of 6 subject simulations") simulate(mgc, cohort; steady_state = opts, replicates = 2, seed = 5, ensemble = EnsembleSerial())
    @test r2 == r2s
    @test count(s -> s.retcode === :SteadyStateFailed, r2.subjects) == 2
    # a cohort of only failures is fine too; the single-subject method still throws
    allbad = [mk("b1", 0.1), mk("b2", 0.2)]
    rb = @test_logs (:warn, r"2 of 2") simulate(mgc, allbad; steady_state = opts)
    @test all(o -> isnan(o.dv), rb.obs)
    @test_throws SteadyStateError simulate(mgc, allbad[1]; η = :zero, steady_state = opts)
    # malformed regimens still throw for the whole call, whatever the other subjects do
    early = Subject("early"; regimen = Regimen(1.0; cmt = :c, ii = 6.0, ss = true, time = -3.0), covariates = (G = -0.3,), obstimes = [1.0])
    @test_throws ArgumentError simulate(mgc, [cohort[1], early]; steady_state = opts)
    # explosive growth ends in a non-finite state, reported as such
    big = PopulationModel(grow; routes = [Route(:c, grow.B)], θ = (g = 5.0,), ηnames = Symbol[], Ω = zeros(0, 0),
        params = (θ, η, c) -> [grow.g => θ.g], observed = (B = grow.B,), error = (B = AdditiveError(1.0),),
        tspan = (0.0, 48.0), solver = Tsit5())
    e2 = with_logger(NullLogger()) do          # the solver warns about the blow-up
        try
            simulate(big, sg; η = :zero); nothing
        catch e
            e
        end
    end
    @test e2 isa SteadyStateError && occursin("no steady state reached after", e2.msg)
    # a stable model that converges too slowly for the limit also fails, with the limit in the message
    slow = pkmodel()
    s = Subject("s"; regimen = Regimen(100.0; cmt = :central, ii = 0.5, ss = true), obstimes = [0.0])
    e3 = try
        simulate(slow, s; η = :zero, steady_state = SteadyStateOptions(maxcycles = 3)); nothing
    catch e
        e
    end
    @test e3 isa SteadyStateError && occursin("after 3 cycles", e3.msg)
    # ... and with enough cycles it converges to the closed form
    r = simulate(slow, s; η = :zero, residual = false, steady_state = SteadyStateOptions(maxcycles = 2000))
    @test rows(r, :A)[1] ≈ 100 / (1 - exp(-K * 0.5)) rtol = 1e-5
    # a steady-state dose before the start of tspan is an error; one after the end is ignored
    @test_throws ArgumentError simulate(slow, Subject("e"; regimen = Regimen(1.0; cmt = :central, ii = 6.0, ss = true, time = -3.0), obstimes = [1.0]); η = :zero)
    late = Subject("l"; regimen = Regimen(100.0; cmt = :central, ii = 6.0, ss = true, time = 500.0), obstimes = [1.0])
    @test rows(simulate(slow, late; η = :zero, residual = false, tspan = (0.0, 100.0)), :A) == [0.0]
    @test_throws ArgumentError SteadyStateOptions(maxcycles = 0)
    @test_throws ArgumentError SteadyStateOptions(rtol = 0.0, atol = 0.0)
end

@testset "steady state in a population: per-subject parameters" begin
    tau, amt = 8.0, 150.0
    m = pkmodel()
    obs_s = [0.0, 2.0, 5.0, 7.9]
    n = 6
    pop = [Subject("p$i"; regimen = Regimen(amt; cmt = :central, ii = tau, ss = true, addl = 1), obstimes = obs_s) for i in 1:n]
    res = simulate(m, pop; seed = 3, residual = false)
    ηs = [r.η.CL for r in res.subjects]
    @test length(unique(ηs)) == n                                  # genuinely different subjects
    peaks = Float64[]
    for (i, s) in enumerate(pop)
        k = 2.0 * exp(ηs[i]) / 10.0
        got = [o.ipred for o in res.obs if o.id == s.id && o.endpoint === :A]
        want = [amt * exp(-k * x) / (1 - exp(-k * tau)) for x in obs_s]
        @test got ≈ want rtol = 1e-6
        push!(peaks, got[1])
    end
    @test length(unique(round.(peaks; digits = 3))) == n          # different parameters, different steady states
    # zero random effects: the deterministic steady state
    z = simulate(m, pop[1:1]; η = :zero, residual = false)
    @test [o.ipred for o in z.obs if o.endpoint === :A] ≈ [amt * exp(-K * x) / (1 - exp(-K * tau)) for x in obs_s] rtol = 1e-6
    @test z.obs == zero_re_prediction(m, pop[1]).obs
    # threaded and serial agree exactly; so do repeated runs; and subject order does not matter
    serial = simulate(m, pop; seed = 3, residual = true, ensemble = EnsembleSerial())
    threaded = simulate(m, pop; seed = 3, residual = true, ensemble = EnsembleThreads())
    @test serial.obs == threaded.obs && serial.subjects == threaded.subjects
    @test Set(simulate(m, reverse(pop); seed = 3, residual = true).obs) == Set(serial.obs)
    # a cohort with and without steady-state subjects, replicates
    mixed = [pop[1], Subject("plain"; regimen = Regimen(amt; cmt = :central, ii = tau, addl = 1), obstimes = obs_s), pop[2]]
    rm = simulate(m, mixed; seed = 4, replicates = 2, residual = false)
    @test length(rm.obs) == 2 * 3 * length(obs_s) * 2                  # 2 reps x 3 subjects x 4 times x 2 endpoints
    plain = [o.ipred for o in rm.obs if o.id == "plain" && o.replicate == 1 && o.endpoint === :A]
    k = 2.0 * exp(rm.subjects[2].η.CL) / 10.0
    @test plain ≈ [amt * exp(-k * x) for x in obs_s] rtol = 1e-6       # starts from zero, no steady state
    # an oral steady-state cohort: depot and central
    oral = [Subject("o$i"; regimen = Regimen(amt; cmt = :depot, ii = tau, ss = true), obstimes = obs_s, control = DoseControl(lag = Dict(:depot => 1.5))) for i in 1:3]
    ro = simulate(m, oral; seed = 9, residual = false)
    for (i, s) in enumerate(oral)
        k = 2.0 * exp(ro.subjects[i].η.CL) / 10.0
        got = [o.ipred for o in ro.obs if o.id == s.id && o.endpoint === :A]
        resp(x) = x < 1.5 ? 0.0 : amt * KA / (KA - k) * (exp(-k * (x - 1.5)) - exp(-KA * (x - 1.5)))
        @test got ≈ [ss_sum(resp, tau, x) for x in obs_s] rtol = 1e-6 atol = 1e-8
    end
end

@testset "steady state uses the covariates at the steady-state time" begin
    tau, amt, T = 6.0, 100.0, 10.0
    m = pkmodel()
    s = Subject("w"; regimen = Regimen(amt; cmt = :central, time = T, ii = tau, ss = true), covariates = (WT = 70.0,),
        timevarying = [4.0 => (WT = 100.0,)], obstimes = [T, T + 2.0])
    r = zero_re_prediction(m, s)
    k = 2.0 * (100.0 / 70)^0.75 / 10.0
    @test rows(r, :A) ≈ [amt / (1 - exp(-k * tau)), amt * exp(-k * 2) / (1 - exp(-k * tau))] rtol = 1e-6
end

@testset "convergence compares componentwise (zero states, atol = 0)" begin
    # the IV bolus never touches the depot: it stays exactly 0
    s = Subject("z"; regimen = Regimen(100.0; cmt = :central, ii = 6.0, ss = true), obstimes = [0.0, 3.0])
    r = zero_re_prediction(pkmodel(), s; steady_state = SteadyStateOptions(atol = 0.0, rtol = 1e-8))
    @test rows(r, :Dep) == [0.0, 0.0]
    @test rows(r, :A) ≈ [100 / (1 - exp(-K * 6)), 100 * exp(-K * 3) / (1 - exp(-K * 6))] rtol = 1e-6
    # atol only is fine as well
    r2 = zero_re_prediction(pkmodel(), s; steady_state = SteadyStateOptions(atol = 1e-9, rtol = 0.0))
    @test rows(r2, :A)[1] ≈ rows(r, :A)[1] rtol = 1e-6
end

@parameters kf = 0.3
@variables Xf(t) = 0.0
@mtkcompile forced = System([D(Xf) ~ -kf * Xf + sin(2π * t / 6.0)], t)
@parameters kd = 0.3
@variables Xd(t) = 0.0 Yd(t)
@mtkcompile decaying = System([D(Xd) ~ -kd * Xd + exp(-0.1 * t), Yd ~ 2 * Xd], t)

@testset "explicit time dependence is refused (unless asserted periodic)" begin
    @test !QSPAtlas.explicit_time_dependence(pk)
    @test QSPAtlas.explicit_time_dependence(forced)
    @test QSPAtlas.explicit_time_dependence(decaying)
    fm = simplemodel(forced, [Route(:c, forced.Xf)], [forced.kf => 0.3]; observed = (X = forced.Xf,))
    dm = simplemodel(decaying, [Route(:c, decaying.Xd)], [decaying.kd => 0.3]; observed = (X = decaying.Xd,))
    reg = Regimen(10.0; cmt = :c, ii = 6.0, ss = true)
    s = Subject("f"; regimen = reg, obstimes = [0.0, 3.0])
    err = try simulate(fm, s; η = :zero); nothing catch e e end
    @test err isa ArgumentError && occursin("explicitly on time", err.msg)
    @test_throws ArgumentError simulate(dm, s; η = :zero)
    @test_throws ArgumentError simulate(fm, [s, Subject("g"; regimen = reg, obstimes = [1.0])])     # a cohort too
    # no steady-state dose: time dependence is fine
    @test simulate(fm, Subject("n"; regimen = Regimen(10.0; cmt = :c), obstimes = [1.0]); η = :zero) isa SimResult
    # asserting periodicity: forcing sin(2 pi t / 6) has period ii = 6 and T = 0; compare with a long run of ordinary doses
    opts = SteadyStateOptions(allow_time_dependent = true)
    r = zero_re_prediction(fm, s; steady_state = opts)
    brute = Subject("b"; regimen = Regimen(10.0; cmt = :c, time = -600.0, ii = 6.0, addl = 100), obstimes = [0.0, 3.0])
    rb = zero_re_prediction(fm, brute; tspan = (-600.0, 12.0))
    @test rows(r, :X) ≈ rows(rb, :X) rtol = 1e-6
end

@testset "accuracy with the solver's DEFAULT tolerances" begin
    m = pkmodel(solve_kwargs = (;))
    tau, amt = 6.0, 100.0
    check_ss(; amt, cmt = :central, tau, model = m, rtol = 1e-5, atol = 1e-4, resp_c = x -> bolus_c(x; amt), label = "default tol bolus")
    check_ss(; amt, cmt = :depot, tau = 12.0, F = 0.8, model = m, rtol = 1e-5, atol = 1e-4, resp_c = x -> oral_c(x; amt = 0.8 * amt),
        resp_d = x -> oral_d(x; amt = 0.8 * amt), label = "default tol oral")
    check_ss(; amt, cmt = :central, tau, duration = 10.0, addl = 2, model = m, rtol = 1e-5, atol = 1e-4, resp_c = x -> infusion_c(x; amt, dur = 10.0), label = "default tol infusion")
    check_ss(; amt, cmt = :depot, tau, lag = 7.0, addl = 1, model = m, rtol = 1e-5, atol = 1e-4, resp_c = x -> oral_c(x; amt, lag = 7.0), label = "default tol lag")
end

@parameters kel = 0.3 kon = 0.05 koff = 0.1
@variables Dr(t) = 0.0 Rc(t) = 10.0 DRc(t) = 0.0
@mtkcompile bind = System([D(Dr) ~ -kel * Dr - kon * Dr * Rc + koff * DRc, D(Rc) ~ -kon * Dr * Rc + koff * DRc,
    D(DRc) ~ kon * Dr * Rc - koff * DRc], t)

@testset "endogenous / conserved state: the initial conditions matter" begin
    # free receptor Rc (initial 10) binds the dosed drug; Rc + DRc = 10 is conserved, so the steady
    # state depends on the initial conditions and cannot be found starting from zero
    m = simplemodel(bind, [Route(:drug, bind.Dr)], [bind.kel => 0.3, bind.kon => 0.05, bind.koff => 0.1];
        observed = (Dr = bind.Dr, R = bind.Rc, DR = bind.DRc))
    s = Subject("b"; regimen = Regimen(20.0; cmt = :drug, ii = 6.0, ss = true), obstimes = [0.0, 3.0, 5.5])
    r = zero_re_prediction(m, s)
    brute = Subject("l"; regimen = Regimen(20.0; cmt = :drug, time = -1800.0, ii = 6.0, addl = 300), obstimes = [0.0, 3.0, 5.5])
    rb = zero_re_prediction(m, brute; tspan = (-1800.0, 12.0))
    for ep in (:Dr, :R, :DR)
        @test rows(r, ep) ≈ rows(rb, ep) rtol = 1e-6
    end
    @test rows(r, :R) .+ rows(r, :DR) ≈ [10.0, 10.0, 10.0] rtol = 1e-8       # conserved total from the ICs
    @test all(>(0), rows(r, :DR))
    # the steady state follows a changed initial condition (total receptor 4)
    m4 = PopulationModel(bind; routes = [Route(:drug, bind.Dr)], θ = (;), ηnames = Symbol[], Ω = zeros(0, 0),
        params = (θ, η, c) -> [bind.kel => 0.3, bind.kon => 0.05, bind.koff => 0.1], init = (θ, η, c) -> [bind.Rc => 4.0],
        observed = (R = bind.Rc, DR = bind.DRc), error = (R = AdditiveError(0.1), DR = AdditiveError(0.1)),
        tspan = (0.0, 100.0), solver = Tsit5(), solve_kwargs = TOL)
    r4 = zero_re_prediction(m4, s)
    @test rows(r4, :R) .+ rows(r4, :DR) ≈ [4.0, 4.0, 4.0] rtol = 1e-8
end

@parameters k10 = 0.2 k12 = 0.15 k21 = 0.1
@variables A1(t) = 0.0 A2(t) = 0.0
@mtkcompile twoc = System([D(A1) ~ -(k10 + k12) * A1 + k21 * A2, D(A2) ~ k12 * A1 - k21 * A2], t)

@testset "two-compartment steady state: both states" begin
    m = simplemodel(twoc, [Route(:c, twoc.A1)], [twoc.k10 => 0.2, twoc.k12 => 0.15, twoc.k21 => 0.1]; observed = (A1 = twoc.A1, A2 = twoc.A2))
    M = [-(0.2 + 0.15) 0.1; 0.15 -0.1]
    tau = 6.0
    b = [100.0, 0.0]
    xpost = (I - exp(M * tau)) \ b                                  # state right after the dose at steady state
    obs = [0.0, 1.0, 3.0, 5.5]
    r = zero_re_prediction(m, Subject("t"; regimen = Regimen(100.0; cmt = :c, ii = tau, ss = true), obstimes = obs))
    want = [exp(M * x) * xpost for x in obs]
    @test rows(r, :A1) ≈ [w[1] for w in want] rtol = 1e-6
    @test rows(r, :A2) ≈ [w[2] for w in want] rtol = 1e-6
    @test all(>(1), rows(r, :A2))                                      # the peripheral state is genuinely loaded
end

@testset "lag and infusion duration exactly integer multiples of ii" begin
    tau, amt = 6.0, 100.0
    check_ss(; amt, cmt = :central, tau, lag = 12.0, addl = 1, resp_c = x -> bolus_c(x; amt, lag = 12.0), label = "bolus lag 2ii")
    check_ss(; amt, cmt = :central, tau, lag = 6.0, addl = 1, resp_c = x -> bolus_c(x; amt, lag = 6.0), label = "bolus lag ii")
    check_ss(; amt, cmt = :depot, tau, lag = 6.0, addl = 1, resp_c = x -> oral_c(x; amt, lag = 6.0), resp_d = x -> oral_d(x; amt, lag = 6.0), label = "oral lag ii")
    check_ss(; amt, cmt = :central, tau, duration = 6.0, addl = 2, resp_c = x -> infusion_c(x; amt, dur = 6.0), label = "infusion duration ii")
    check_ss(; amt, cmt = :central, tau, duration = 12.0, addl = 2, resp_c = x -> infusion_c(x; amt, dur = 12.0), label = "infusion duration 2ii")
    check_ss(; amt, cmt = :central, tau, lag = 6.0, duration = 6.0, addl = 2, resp_c = x -> infusion_c(x; amt, dur = 6.0, lag = 6.0), label = "infusion lag ii, duration ii")
    check_ss(; amt, cmt = :central, tau, lag = 12.0, duration = 12.0, addl = 1, T = 6.0, resp_c = x -> infusion_c(x; amt, dur = 12.0, lag = 12.0), label = "infusion lag 2ii, duration 2ii")
end

@parameters kb = 0.2 Rb = 0.5
@variables Xb(t) = 0.0
@mtkcompile bsys = System([D(Xb) ~ -kb * Xb + Rb], t)

@testset "steady-state infusion on top of a nonzero rate baseline" begin
    routes = [Route(:c, bsys.Xb; rate = bsys.Rb)]
    m = simplemodel(bsys, routes, [bsys.kb => 0.2]; observed = (X = bsys.Xb,))
    amt, tau = 100.0, 6.0
    base = 0.5 / 0.2                                                   # level from the background input alone
    for (dur, addl) in ((2.0, 1), (10.0, 2))
        s = Subject("n"; regimen = Regimen(amt; cmt = :c, ii = tau, ss = true, duration = dur, addl), obstimes = OBS_S)
        r = zero_re_prediction(m, s)
        want = [base + ss_sum(x -> infusion_c(x; amt, dur), tau, x; addl) for x in OBS_S]
        @test rows(r, :X) ≈ want rtol = 1e-6
    end
    # no dose at all after the last infusion: the background level is the limit
    s = Subject("n2"; regimen = Regimen(amt; cmt = :c, ii = tau, ss = true, duration = 2.0), obstimes = [0.0, 200.0])
    @test rows(zero_re_prediction(m, s), :X)[2] ≈ base rtol = 1e-6
end

@testset "init -> reinit! keeps the rate baseline (guards the SS cycle loop)" begin
    routes = [Route(:c, bsys.Xb; rate = bsys.Rb)]
    ev = resolve(Regimen(20.0; cmt = :c, time = 1.0, duration = 10.0), routes)      # extra rate 2 on [1, 11)
    cb, ts = dosing_callback(bsys, ev, routes)
    prob = ODEProblem(bsys, [], (0.0, 20.0))
    step(x0, r, dt) = r / 0.2 + (x0 - r / 0.2) * exp(-0.2 * dt)
    x1 = step(0.0, 0.5, 1.0); x11 = step(x1, 2.5, 10.0); x20 = step(x11, 0.5, 9.0)
    integ = init(prob, Tsit5(); callback = cb, tstops = ts, abstol = 1e-12, reltol = 1e-12, save_on = false)
    for round in 1:3
        # the first initialization already switched the infusion state on; reinit! must not take that as the baseline
        reinit!(integ, [0.0]; t0 = 0.0, tf = 20.0)
        solve!(integ)
        @test integ.u[1] ≈ x20 rtol = 1e-7
        @test integ.ps[bsys.Rb] == 0.5
    end
    reinit!(integ)                                                    # plain reinit! as well
    solve!(integ)
    @test integ.u[1] ≈ x20 rtol = 1e-7
    @test integ.ps[bsys.Rb] == 0.5
end

end
