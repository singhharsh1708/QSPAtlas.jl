# NONMEM-style datasets (read_nmdata / write_nmdata), record order (:pre / :post observations),
# NOCB and LOCF covariates, EVID 2, SS=2, constant-infusion steady state and modelled rates and
# durations. Own module, so helper names cannot clash with other test files. Expected values are
# closed forms of the one-compartment model (amount A, A' = -k A + input, k = CL / V).
module DatasetTests

using QSPAtlas, Test
using QSPAtlas: covariates_at, with_obstimes
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

const TOL = (abstol = 1e-11, reltol = 1e-11)

@parameters CL = 2.0 V = 10.0 Rin = 0.0 ka = 1.1 KIN = 10.0 KOUT = 0.1
@variables A(t) = 0.0 Ad(t) = 0.0 R(t) = 100.0
@mtkcompile pk = System([D(Ad) ~ -ka * Ad, D(A) ~ ka * Ad - CL / V * A + Rin, D(R) ~ KIN - KOUT * R], t)
const ROUTES = [Route(:depot, pk.Ad), Route(:central, pk.A; rate = pk.Rin)]
const CMT = Dict(1 => :depot, 2 => :central)

# CL = 2 (WT / 70): k = 0.2 at WT 70
model(; tspan = (0.0, 100.0)) = PopulationModel(pk; routes = ROUTES, θ = (CL = 2.0, V = 10.0), ηnames = Symbol[], Ω = zeros(0, 0),
    params = (θ, η, c) -> [pk.CL => θ.CL * get(c, :WT, 70.0) / 70, pk.V => θ.V],
    observed = (A = pk.A, C = pk.A / pk.V, Dep = pk.Ad, R = pk.R),
    error = (A = AdditiveError(0.1), C = AdditiveError(0.1), Dep = AdditiveError(0.1), R = AdditiveError(0.1)),
    tspan, solver = Tsit5(), solve_kwargs = TOL)
const M = model()
const K = 0.2

# one record as a NamedTuple with every column ("." = missing)
rec(id, time, evid; amt = ".", cmt = 2, rate = ".", ss = ".", ii = ".", addl = ".", mdv = evid == 0 ? 0 : 1, dv = ".", wt = 70) =
    (ID = id, TIME = time, EVID = evid, AMT = amt, CMT = cmt, RATE = rate, SS = ss, II = ii, ADDL = addl, MDV = mdv, DV = dv, WT = wt)
obs(id, time; kw...) = rec(id, time, 0; kw...)
dose(id, time, amt; kw...) = rec(id, time, 1; amt, kw...)

pred(res, ep = :A) = [o.ipred for o in res.obs if o.endpoint === ep]
sim(data; kw...) = simulate(M, data.subjects; η = :zero, residual = false, kw...)

@testset "record order: :pre and :post observations at a dose time" begin
    rows = [obs(1, 0.0), dose(1, 0.0, 100; ii = 12, addl = 1), obs(1, 0.0), obs(1, 6.0), obs(1, 12.0),
            obs(1, 24.0), dose(1, 24.0, 100), obs(1, 24.0), obs(1, 30.0)]
    data = read_nmdata(rows; cmt = CMT, endpoints = :all)
    @test [o.order for o in data.obs] == [:pre, :post, :post, :post, :pre, :post, :post]
    a(t) = 100exp(-K * t) + (t >= 12 ? 100exp(-K * (t - 12)) : 0.0)
    want = [0.0, 100.0, a(6.0), a(12.0), a(24.0), a(24.0) + 100, a(30.0) + 100exp(-K * 6)]
    @test pred(sim(data)) ≈ want rtol = 1e-8 atol = 1e-9
    # the ADDL dose at 12 comes before the observation at 12 (post), as in mrgsolve's default
    @test pred(sim(data))[4] ≈ 100exp(-K * 12) + 100 rtol = 1e-8

    # the same through the Subject API: default :post keeps the old behaviour
    reg = Regimen(100.0; cmt = :central, ii = 12.0, addl = 1)
    post = zero_re_prediction(M, Subject("p"; regimen = reg, obstimes = [12.0]))
    tagged = zero_re_prediction(M, Subject("q"; regimen = reg, obstimes = [12.0, 12.0], obsorder = [:post, :pre]))
    @test pred(post) ≈ [100exp(-K * 12) + 100] rtol = 1e-8
    @test pred(tagged) ≈ [100exp(-K * 12), 100exp(-K * 12) + 100] rtol = 1e-8       # stored :pre first
    # pre at the initial time and at a reset
    s = Subject("r"; regimen = Regimen(Regimen(100.0; cmt = :central), Reset(5.0)), obstimes = [0.0, 0.0, 5.0, 5.0],
        obsorder = [:pre, :post, :pre, :post])
    @test pred(zero_re_prediction(M, s)) ≈ [0.0, 100.0, 100exp(-K * 5), 0.0] atol = 1e-8
    @test_throws ArgumentError Subject("x"; regimen = reg, obstimes = [1.0], obsorder = [:before])
    @test_throws ArgumentError Subject("x"; regimen = reg, obstimes = [1.0, 2.0], obsorder = [:pre])
end

@testset "endpoints per observation" begin
    rows = [dose(1, 0.0, 100), obs(1, 1.0; cmt = 2), obs(1, 2.0; cmt = 1), obs(1, 3.0; cmt = 0)]
    data = read_nmdata(rows; cmt = CMT, endpoints = Dict(2 => :C, 1 => :Dep))
    @test [o.endpoint for o in data.obs] == [:C, :Dep, nothing]
    res = sim(data)
    @test [(o.time, o.endpoint) for o in res.obs] == [(3.0, :A), (1.0, :C), (3.0, :C), (2.0, :Dep), (3.0, :Dep), (3.0, :R)]
    @test_throws ArgumentError read_nmdata(rows; cmt = CMT)                      # numeric obs CMT needs a mapping (or :all)
    bad = Subject("b"; regimen = Regimen(1.0; cmt = :central), obstimes = [1.0], endpoints = [:nope])
    @test_throws ArgumentError simulate(M, bad; η = :zero)
end

@testset "round trip read / write" begin
    rows = [dose(1, 0.0, 100; ii = 12, addl = 2, wt = 60), obs(1, 0.0; wt = 60, dv = 9.5), obs(1, 6.0; wt = 80, dv = 3.25),
            rec(1, 9.0, 2; cmt = ".", wt = 90), obs(1, 12.0; wt = 90), rec(1, 13.0, 3; wt = 90), dose(1, 14.0, 50; cmt = 2, rate = 25, wt = 90),
            obs(1, 15.0; wt = 90), dose(2, 0.0, 0; rate = 5, ss = 1), dose(2, 0.0, 30; ss = 2, ii = 24), obs(2, 0.0), obs(2, 4.0),
            dose(3, 0.0, 100; rate = -2), obs(3, 1.0), rec(3, 2.0, 4; amt = 10, cmt = "central"), obs(3, 2.0)]
    ctl = DoseControl(duration = Dict(:central => 4.0))
    d1 = read_nmdata(rows; cmt = CMT, endpoints = :all, control = ctl)
    dir = mktempdir()
    f1 = write_nmdata(joinpath(dir, "a.csv"), d1)
    d2 = read_nmdata(f1; cmt = CMT, endpoints = :all, control = ctl)
    @test d2.columns == d1.columns
    @test all(isequal.(d2.values, d1.values))
    @test isequal(d2.obs, d1.obs)
    for (a, b) in zip(d1.subjects, d2.subjects)
        @test a.id == b.id && a.obstimes == b.obstimes && a.obsorder == b.obsorder && a.obsendpoints == b.obsendpoints
        @test a.covariates == b.covariates && a.covariate_times == b.covariate_times && a.covariate_values == b.covariate_values
        @test doses(a.regimen) == doses(b.regimen) && reset_times(a.regimen) == reset_times(b.regimen)
    end
    s1 = d1.subjects[1]
    @test s1.covariate_rule === :nocb && s1.covariates.WT == 60
    @test s1.covariate_times == [0.0, 6.0, 15.0] && [v.WT for v in s1.covariate_values] == [60.0, 80.0, 90.0]
    @test reset_times(s1.regimen) == [13.0]
    @test [d.ss for d in doses(d1.subjects[2].regimen)] == [1, 2]
    @test doses(d1.subjects[3].regimen)[1].fixed === :model_duration
    @test reset_times(d1.subjects[3].regimen) == [2.0]
    @test isequal([o.dv for o in d1.obs[1:2]], [9.5, 3.25])
    # simulated DV written back record by record
    res = simulate(M, d1.subjects; η = :zero, residual = false)
    f2 = write_nmdata(joinpath(dir, "b.csv"), d1; result = res, value = :ipred, endpoint = :C)
    d3 = read_nmdata(f2; cmt = CMT, endpoints = :all, control = ctl)
    want = Float64[]
    for s in d1.subjects
        append!(want, [o.ipred for o in res.obs if o.id == s.id && o.endpoint === :C])
    end
    @test [o.dv for o in d3.obs] == want
    @test_throws ArgumentError write_nmdata(joinpath(dir, "c.csv"), d1; result = res)     # which endpoint?
    # a whitespace-separated file with lower-case names reads the same
    f3 = joinpath(dir, "ws.txt")
    write(f3, "id time evid amt cmt\n1 0 1 100 2\n1 1 0 . 2\n")
    d4 = read_nmdata(f3; cmt = CMT, endpoints = :all)
    @test pred(sim(d4)) ≈ [100exp(-K)] rtol = 1e-8
end

@testset "covariates: NOCB vs LOCF" begin
    # CL = 2 WT / 70: k = 0.2 WT / 70
    rows = [dose(1, 0.0, 100; wt = 70), obs(1, 4.0; wt = 140), obs(1, 10.0; wt = 35), obs(1, 20.0; wt = 35)]
    k(w) = 0.2w / 70
    nocb = read_nmdata(rows; cmt = CMT, endpoints = :all)
    locf = read_nmdata(rows; cmt = CMT, endpoints = :all, covariate_rule = :locf)
    # NOCB: (0, 4] at WT 140, (4, 10] and (10, 20] at 35
    @test pred(sim(nocb)) ≈ [100exp(-k(140) * 4), 100exp(-k(140) * 4 - k(35) * 6), 100exp(-k(140) * 4 - k(35) * 16)] rtol = 1e-8
    # LOCF: [0, 4) at 70, [4, 10) at 140, then 35
    @test pred(sim(locf)) ≈ [100exp(-k(70) * 4), 100exp(-k(70) * 4 - k(140) * 6), 100exp(-k(70) * 4 - k(140) * 6 - k(35) * 10)] rtol = 1e-8
    # the observable at a record uses that record's value under both rules (C = A / V does not
    # depend on WT here; covariates_at shows the value in force)
    s = nocb.subjects[1]
    # the dose record at 0 (WT 70) ends the first interval: 70 up to 0, 140 on (0, 4]
    @test [covariates_at(s, x).WT for x in (-1.0, 0.0, 2.0, 4.0, 5.0, 10.0, 30.0)] == [70, 70, 140, 140, 35, 35, 35]
    sl = locf.subjects[1]
    @test [covariates_at(sl, x).WT for x in (0.0, 2.0, 4.0, 5.0, 10.0)] == [70, 70, 140, 140, 35]
    # the Subject API with covariate_rule = :nocb gives the same
    sub = Subject("n"; regimen = Regimen(100.0; cmt = :central), obstimes = [4.0, 10.0], covariates = (WT = 140.0,),
        timevarying = [4.0 => (WT = 140.0,), 10.0 => (WT = 35.0,)], covariate_rule = :nocb)
    @test pred(zero_re_prediction(M, sub)) ≈ pred(sim(nocb))[1:2] rtol = 1e-10
    @test with_obstimes(sub, [1.0]).covariate_rule === :nocb
    @test_throws ArgumentError Subject("n"; regimen = Regimen(1.0; cmt = :central), obstimes = [1.0], covariate_rule = :nobc)
end

@testset "EVID 2: other-type records" begin
    k(w) = 0.2w / 70
    rows = [dose(1, 0.0, 100), obs(1, 2.0), rec(1, 6.0, 2; wt = 140, cmt = "."), obs(1, 8.0; wt = 140),
            rec(1, 10.0, 2; wt = 35, cmt = 2), obs(1, 14.0; wt = 35)]
    data = read_nmdata(rows; cmt = CMT, endpoints = :all)
    @test length(data.obs) == 3                                   # no output for EVID 2 records
    # NOCB: (2, 6] and (6, 8] at 140 (the EVID 2 record carries the change), (8, 10] and (10, 14] at 35
    a8 = 100exp(-k(70) * 2 - k(140) * 6)
    @test pred(sim(data)) ≈ [100exp(-k(70) * 2), a8, a8 * exp(-k(35) * 6)] rtol = 1e-8
    # EVID 2 is derived from MDV = 1 and AMT = 0 when there is no EVID column
    nod = [(ID = 1, TIME = 0.0, AMT = 100, CMT = 2, MDV = 1, WT = 70), (ID = 1, TIME = 2.0, AMT = 0, CMT = 2, MDV = 0, WT = 70),
           (ID = 1, TIME = 6.0, AMT = 0, CMT = 2, MDV = 1, WT = 140), (ID = 1, TIME = 8.0, AMT = 0, CMT = 2, MDV = 0, WT = 140)]
    @test pred(sim(read_nmdata(nod; cmt = CMT, endpoints = :all))) ≈ pred(sim(data))[1:2] rtol = 1e-10
    @test_throws ArgumentError read_nmdata([dose(1, 0.0, 100), rec(1, 1.0, 2; mdv = 0)]; cmt = CMT, endpoints = :all)
    @test_throws ArgumentError read_nmdata([dose(1, 0.0, 100), rec(1, 1.0, 2; amt = 5)]; cmt = CMT, endpoints = :all)
    @test_throws ArgumentError read_nmdata([dose(1, 0.0, 100), rec(1, 1.0, 2; cmt = -2)]; cmt = CMT, endpoints = :all)
end

@testset "SS=2 against explicit doses" begin
    # NONMEM Users Guide example: 10 every morning (8 h), 15 every evening (20 h), to steady state
    rows = [dose(1, 8.0, 10; ss = 1, ii = 24), obs(1, 14.0), obs(1, 20.0), dose(1, 20.0, 15; ss = 2, ii = 24), obs(1, 20.0),
            obs(1, 26.0), obs(1, 32.0), obs(1, 44.0)]
    data = read_nmdata(rows; cmt = CMT, endpoints = :all)
    got = pred(sim(data))
    # explicit doses from long before: 10 at 8 - 24 j, 15 at 20 - 24 j (j >= 0)
    a(t; before20 = false) = sum(10exp(-K * (t - (8 - 24j))) for j in 0:200 if 8 - 24j <= t) +
        sum(15exp(-K * (t - (20 - 24j))) for j in (before20 ? 1 : 0):200 if 20 - 24j <= t)
    # before the SS=2 record only the SS=1 regimen is there (NONMEM: "will reflect only the first steady-state dosage")
    first(t) = sum(10exp(-K * (t - (8 - 24j))) for j in 0:200 if 8 - 24j <= t)
    @test got ≈ [first(14.0), first(20.0), a(20.0), a(26.0), a(32.0), a(44.0)] rtol = 1e-8
    # second form of the same example: 10 every 12 h plus 5 every 24 h, both at 20
    rows2 = [dose(1, 20.0, 10; ss = 1, ii = 12), dose(1, 20.0, 5; ss = 2, ii = 24), obs(1, 20.0), obs(1, 26.0), obs(1, 32.0), obs(1, 44.0)]
    @test pred(sim(read_nmdata(rows2; cmt = CMT, endpoints = :all))) ≈ got[3:end] rtol = 1e-8
    # an oral SS=2 dose added to a non-steady-state history; the endogenous R is not counted twice
    reg = Regimen(Regimen(100.0; cmt = :central, time = 0.0), Regimen(50.0; cmt = :depot, time = 10.0, ii = 24.0, ss = 2))
    res = zero_re_prediction(M, Subject("o"; regimen = reg, obstimes = [10.0, 15.0]))
    expl = Regimen(Regimen(100.0; cmt = :central, time = 0.0), Regimen(50.0; cmt = :depot, time = 10.0 - 24 * 60, ii = 24.0, addl = 60))
    ref = zero_re_prediction(model(tspan = (-1500.0, 100.0)), Subject("e"; regimen = expl, obstimes = [10.0, 15.0]))
    @test pred(res) ≈ pred(ref) rtol = 1e-7
    @test pred(res, :R) ≈ [100.0, 100.0] rtol = 1e-10
    @test_throws ArgumentError Regimen(1.0; cmt = :central, ss = 3, ii = 12.0)
end

@testset "constant infusion at steady state" begin
    # AMT 0, RATE 5, SS 1: amount R / k, concentration R / CL; the infusion ends at the record (NONMEM)
    rows = [dose(1, 0.0, 0; rate = 5, ss = 1), obs(1, 0.0), obs(1, 3.0),
            dose(2, 0.0, 0; rate = 5, ss = 1), dose(2, 0.0, 500; rate = 5), obs(2, 0.0), obs(2, 3.0), obs(2, 50.0)]
    data = read_nmdata(rows; cmt = CMT, endpoints = :all)
    r = sim(data)
    @test pred(r, :C)[1] ≈ 5 / 2.0 rtol = 1e-8                               # Css = R / CL
    @test pred(r)[2] ≈ 25exp(-K * 3) rtol = 1e-8
    @test pred(r, :C)[3:5] ≈ fill(2.5, 3) rtol = 1e-8                          # kept running by an explicit infusion
    # at a later time, after other doses (SS=1 replaces them), with WT 140 (CL 4) and modelled rate
    reg = Regimen(Regimen(100.0; cmt = :central), Regimen(0.0; time = 7.0, cmt = :central, rate = :model, ss = true))
    s = Subject("m"; regimen = reg, obstimes = [7.0], covariates = (WT = 140.0,), control = DoseControl(rate = Dict(:central => 8.0)))
    @test pred(zero_re_prediction(M, s), :C) ≈ [8 / 4.0] rtol = 1e-8
    # bioavailability does not apply (NONMEM), lag has no effect
    s2 = Subject("f"; regimen = Regimen(0.0; cmt = :central, rate = 5.0, ss = 1), obstimes = [0.0],
        control = DoseControl(bioav = Dict(:central => 0.5), lag = Dict(:central => 3.0)))
    @test pred(zero_re_prediction(M, s2), :C) ≈ [2.5] rtol = 1e-8
    # SS=2 constant infusion adds its steady state
    s3 = Subject("a"; regimen = Regimen(Regimen(100.0; cmt = :central, time = 0.0), Regimen(0.0; time = 5.0, cmt = :central, rate = 5.0, ss = 2)),
        obstimes = [5.0])
    @test pred(zero_re_prediction(M, s3)) ≈ [100exp(-K * 5) + 25] rtol = 1e-8
    @test_throws ArgumentError Regimen(0.0; cmt = :central, rate = 5.0, ss = true, ii = 12.0)
    @test_throws ArgumentError Regimen(0.0; cmt = :central, ss = true, ii = 12.0)
    @test_throws ArgumentError read_nmdata([dose(1, 0.0, 0; rate = 5, ss = 1, ii = 12)]; cmt = CMT)
    @test_throws ArgumentError read_nmdata([dose(1, 0.0, 0; rate = 5)]; cmt = CMT)       # AMT 0 without SS
end

@testset "EVID inferred without an EVID column (AMT, RATE, SS)" begin
    noevid(r) = Base.structdiff(r, NamedTuple{(:EVID,)})
    nomdv(r) = Base.structdiff(r, NamedTuple{(:EVID, :MDV)})
    # a steady-state infusion record (AMT 0, RATE 5, SS 1) is a dose; observations follow it
    rows = [dose(1, 0.0, 0; rate = 5, ss = 1), obs(1, 0.0), obs(1, 3.0)]
    want = pred(sim(read_nmdata(rows; cmt = CMT, endpoints = :all)))
    @test want ≈ [25.0, 25exp(-K * 3)] rtol = 1e-8                     # amount R / k, then the infusion ends
    for f in (noevid, nomdv)
        d = read_nmdata(f.(rows); cmt = CMT, endpoints = :all)
        @test length(d.obs) == 2
        @test pred(sim(d)) ≈ want rtol = 1e-12
    end
    # modelled rate (RATE -1) with AMT 0 and SS 1, and an SS=2 infusion record: doses as well
    ctl = DoseControl(rate = Dict(:central => 5.0))
    rm = [dose(1, 0.0, 0; rate = -1, ss = 1), obs(1, 3.0)]
    @test pred(sim(read_nmdata(nomdv.(rm); cmt = CMT, endpoints = :all, control = ctl))) ≈ [25exp(-K * 3)] rtol = 1e-8
    r2 = [dose(1, 0.0, 100), obs(1, 5.0), dose(1, 5.0, 0; rate = 5, ss = 2), obs(1, 5.0)]
    @test pred(sim(read_nmdata(noevid.(r2); cmt = CMT, endpoints = :all))) ≈ [100exp(-K * 5), 100exp(-K * 5) + 25] rtol = 1e-8
    # MDV 1 and AMT, RATE, SS all 0: still an other-type record (no output row)
    @test length(read_nmdata(noevid.([dose(1, 0.0, 100), rec(1, 2.0, 2), obs(1, 3.0)]); cmt = CMT, endpoints = :all).obs) == 1
    # a dose record with MDV 0 is a dose, and MDV must be 1 on it
    @test_throws ArgumentError read_nmdata(noevid.([dose(1, 0.0, 0; rate = 5, ss = 1, mdv = 0)]); cmt = CMT)
    # RATE or SS alone makes a dose record, which then needs a valid dose
    @test_throws ArgumentError read_nmdata(nomdv.([dose(1, 0.0, 0; rate = 5)]); cmt = CMT)
    @test_throws ArgumentError read_nmdata(nomdv.([dose(1, 0.0, 0; ss = 1, ii = 12)]); cmt = CMT)
end

@testset "RATE -1 / -2: modelled rate and duration" begin
    inf(t; R, dur) = t <= dur ? R / K * (1 - exp(-K * t)) : R / K * (1 - exp(-K * dur)) * exp(-K * (t - dur))
    rows = [dose(1, 0.0, 100; rate = -1), obs(1, 1.0), obs(1, 10.0), dose(2, 0.0, 100; rate = -2), obs(2, 1.0), obs(2, 10.0)]
    ctl = DoseControl(bioav = Dict(:central => 0.5), rate = Dict(:central => 20.0), duration = Dict(:central => 5.0))
    r = sim(read_nmdata(rows; cmt = CMT, endpoints = :all, control = ctl))
    # -1: rate 20 kept, F shortens the duration to 2.5; -2: duration 5 kept, F halves the rate to 10
    @test pred(r) ≈ [inf(1.0; R = 20, dur = 2.5), inf(10.0; R = 20, dur = 2.5), inf(1.0; R = 10, dur = 5), inf(10.0; R = 10, dur = 5)] rtol = 1e-8
    # control as a function of id and covariates
    f = (id, c) -> DoseControl(rate = Dict(:central => c.WT / 3.5))
    r2 = sim(read_nmdata([dose(1, 0.0, 100; rate = -1), obs(1, 1.0)]; cmt = CMT, endpoints = :all, control = f))
    @test pred(r2) ≈ [inf(1.0; R = 20, dur = 5)] rtol = 1e-8
    # a modelled value the control does not give: a clear error when simulated
    d = read_nmdata([dose(1, 0.0, 100; rate = -2), obs(1, 1.0)]; cmt = CMT, endpoints = :all)
    @test_throws ArgumentError sim(d)
    @test_throws ArgumentError read_nmdata([dose(1, 0.0, 100; rate = -3)]; cmt = CMT)
    @test_throws ArgumentError Dose(0.0, 1.0; cmt = :a, rate = :model, duration = :model)
    @test_throws ArgumentError Dose(0.0, 1.0; cmt = :a, rate = :fast)
    @test doses(Regimen(100.0; cmt = :central, duration = :model))[1].fixed === :model_duration
end

@testset "unsupported input is an error" begin
    E(rows; kw...) = @test_throws ArgumentError read_nmdata(rows; cmt = CMT, endpoints = :all, kw...)
    E([merge(obs(1, 0.0), (PCMT = 1,))])                                            # unsupported data item
    E([obs(1, 0.0)]; covariates = Symbol[])                                          # WT neither covariate nor dropped
    @test read_nmdata([obs(1, 0.0)]; cmt = CMT, endpoints = :all, drop = [:WT]).covariates == Symbol[]
    E([rec(1, 0.0, 5)])                                                              # EVID 5
    E([dose(1, 0.0, 10; ss = 3, ii = 12)])                                           # SS 3
    E([obs(1, 0.0; mdv = 1)])                                                        # observation with MDV 1
    E([dose(1, 0.0, 10; mdv = 0)])                                                   # dose with MDV 0
    E([obs(1, 2.0), obs(1, 1.0)])                                                    # TIME decreases
    E([obs(1, 0.0), obs(2, 0.0), obs(1, 1.0)])                                       # ID reappears
    E([dose(1, 0.0, 10), obs(1, 0.0), dose(1, 0.0, 10)])                             # observation between two doses at one time
    E([dose(1, 0.0, 10), rec(1, 0.0, 3)])                                            # reset after a dose at the same time
    E([dose(1, 0.0, 10), dose(1, 0.0, 10; ss = 1, ii = 12)])                         # SS=1 after a dose at the same time
    E([dose(1, 0.0, 10; ii = 12, addl = 1), obs(1, 12.0), dose(1, 12.0, 5)])         # pre observation at an ADDL time
    E([obs(1, 0.0; wt = 70), obs(1, 0.0; wt = 80)])                                  # covariates differ at one time
    E([obs(1, 0.0; wt = ".")])                                                       # missing covariate
    E([obs(1, 0.0; wt = "heavy")])                                                   # text covariate
    E([dose(1, 0.0, 10; cmt = -2)])                                                  # turn a compartment off
    E([dose(1, 0.0, 10; cmt = 3)])                                                   # CMT without a route
    E([dose(1, 0.0, 10; cmt = ".")])                                                 # dose without CMT
    E([dose(1, 0.0, 0)])                                                             # AMT 0 bolus
    E([dose(1, 0.0, 10; ii = 12)])                                                   # II without SS or ADDL
    E([dose(1, 0.0, 10; addl = 2)])                                                  # ADDL without II
    E([obs(1, 0.0; amt = 5)])                                                        # AMT on an observation
    E([rec(1, 0.0, 4; amt = 10, ss = 1, ii = 12)])                                   # SS on EVID 4
    E([(ID = 1, EVID = 0)])                                                          # no TIME
    E([(ID = 1, TIME = 0.0, time = 1.0)])                                            # TIME twice
    E([dose(1, 0.0, 10; ss = 1, ii = 12), dose(1, 24.0, 10; ss = 1, ii = 12)])       # two SS=1 times
    dir = mktempdir()
    f = joinpath(dir, "c.csv")
    write(f, "ID,TIME,EVID\n1,0,0\n#1,1,0\n")
    @test_throws ArgumentError read_nmdata(f; endpoints = :all)                      # comment record
    write(f, "ID,TIME,EVID\n1,0,0\n1,1\n")
    @test_throws ArgumentError read_nmdata(f; endpoints = :all)                      # short record
    @test_throws ArgumentError read_nmdata(joinpath(dir, "none.csv"))
    @test_throws ArgumentError read_nmdata([obs(1, 0.0)]; covariate_rule = :last)
end

end # module DatasetTests
