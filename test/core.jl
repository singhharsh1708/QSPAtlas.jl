# Unit tests for the scoring core, using an analytic one-compartment model.
using QSPAtlas, Test
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "reference round trip" begin
    r = ReferenceCurve("c1", "Fig 1", [0.0, 1.0], [1.0, 2.0], Axis("time", "h", :linear, 0, 10),
        Axis("conc", "mg/L", :log, 0.01, 100), :curve, Dict("source" => "test"))
    p = write_reference(tempname() * ".csv", r)
    r2 = load_reference(p)
    @test r2.id == "c1" && r2.x == r.x && r2.y == r.y && r2.yaxis.scale === :log && r2.meta["source"] == "test"
    @test_throws ArgumentError Axis("x", "u", :log, 0, 1)
end

@testset "one-compartment bolus, scored against its analytic solution" begin
    @parameters k = 0.3
    @variables C(t) = 0.0
    @mtkcompile sys = System([D(C) ~ -k * C], t)
    times = repeated_doses(0.0, 6.0, 3)
    cb = bolus_callback(sys, sys.C, times, 10.0)
    prob = ODEProblem(sys, [], (0.0, 24.0); callback = cb)
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    exact(tt) = sum(10.0 * exp(-0.3 * (tt - td)) for td in times if td <= tt; init = 0.0)
    xs = collect(0.5:0.5:23.5)
    ref = ReferenceCurve("conc", "analytic", xs, exact.(xs), Axis("time", "h", :linear, 0, 24), Axis("C", "mg/L", :log, 0.01, 100), :curve, Dict{String, String}())
    grid = collect(range(0, 24; length = 2000))
    good = run_curves(_ -> (grid, sol(grid; idxs = sys.C).u), [ref], "as published")
    @test QSPAtlas.matched(good)
    bad = run_curves(_ -> (grid, 1.3 .* sol(grid; idxs = sys.C).u), [ref], "as published")
    @test !QSPAtlas.matched(bad) && QSPAtlas.ran(bad)
    broken = run_curves(_ -> throw(NotRunnable("Eq 3 uses an undefined symbol")), [ref], "as published")
    @test_throws ErrorException run_curves(_ -> error("a bug in the model code"), [ref], "as published")
    # a published run with gap fills can match, but is never graded A
    filled = run_curves(_ -> (grid, sol(grid; idxs = sys.C).u), [ref], "as published";
        fixes = [Fix("mw", :missing, "molecular weight"; evidence = "IgG", basis = :standard, stage = :before_comparison)])
    @test grade(Reproduction(id = "x", citation = "", doi = "", published = filled)) == 'B'
    @test grade_strict(Reproduction(id = "x", citation = "", doi = "", published = filled)) == 'D'
    @test grade_strict(Reproduction(id = "x", citation = "", doi = "", published = good)) == 'A'
    @test_throws ArgumentError Fix("p0", :from_target, "start value"; evidence = "fig", basis = :paper)
    @test !QSPAtlas.ran(broken)
    @test grade(Reproduction(id = "x", citation = "", doi = "", published = good)) == 'A'
    @test grade(Reproduction(id = "x", citation = "", doi = "", published = bad, fixed = good)) == 'B'
    @test grade(Reproduction(id = "x", citation = "", doi = "", published = bad)) == 'C'
    @test grade(Reproduction(id = "x", citation = "", doi = "", published = broken)) == 'D'
    @test grade(Reproduction(id = "x", citation = "", doi = "", published = good, blocker = "missing")) == 'D'
    rep = Reproduction(id = "x", citation = "c", doi = "d", published = bad, fixed = good)
    s = read_score(write_score(tempname() * ".toml", rep))
    @test s["grade"] == "B" && s["fixed"]["curves"][1]["passed"]
    @test isfile(plot_reproduction(rep, tempname() * ".png"))
    wp = benchmark_problem(BenchmarkSetup(prob; callback = cb, saveat = xs), [Dict(:alg => Tsit5()), Dict(:alg => Rodas5P())];
        reference_alg = Vern9(), abstols = [1e-6, 1e-8], reltols = [1e-4, 1e-6], numruns = 2)
    @test length(wp.wps) == 2
end

@testset "bars are compared vertically at their own position" begin
    ax = Axis("bar", "-", :linear, 0, 4)
    ref = ReferenceCurve("b", "bars", [1.0, 2.0, 3.0], [0.0, -1.0, -2.0], ax, Axis("v", "-", :linear, -3, 1), :bars, Dict{String, String}())
    @test passes(compare_curve(ref, [1.0, 2.0, 3.0], [0.0, -1.0, -2.0]))
    # swapping two bars must fail even though each value appears somewhere
    @test !passes(compare_curve(ref, [1.0, 2.0, 3.0], [0.0, -2.0, -1.0]))
    # a missing or failed bar is a maximum error, not a pass
    c = compare_curve(ref, [1.0, 2.0], [0.0, -1.0])
    @test c.d[3] == 1.0 && !passes(c)
    @test compare_curve(ref, [1.0, 2.0, 3.0], [0.0, NaN, -2.0]).d[2] == 1.0
end

@testset "gaps in a simulation are not bridged" begin
    ax = Axis("t", "-", :linear, 0, 2)
    ref = ReferenceCurve("g", "gap", [1.0], [1.0], ax, Axis("y", "-", :linear, 0, 2), :curve, Dict{String, String}())
    # bridging the gap would give the line y = x through (1, 1); with the gap kept, the
    # nearest simulated points are (0.1, 0.1) and (1.9, 1.9), far from (1, 1)
    c = compare_curve(ref, [0.0, 0.1, 1.0, 1.9, 2.0], [0.0, 0.1, NaN, 1.9, 2.0])
    @test c.d[1] > 0.3 && !passes(c)
    @test_throws ArgumentError compare_curve(ref, [0.0, 1.0, 2.0], [NaN, 1.0, NaN])
end

@testset "figure grade needs curves of that figure" begin
    ax = Axis("t", "-", :linear, 0, 1)
    r1 = ReferenceCurve("fig1_a", "f1", [0.5], [0.5], ax, ax, :curve, Dict{String, String}())
    r2 = ReferenceCurve("fig2_a", "f2", [0.5], [0.5], ax, ax, :curve, Dict{String, String}())
    good = run_curves(_ -> ([0.0, 1.0], [0.0, 1.0]), [r1], "as published")          # only fig1
    both = run_curves(_ -> ([0.0, 1.0], [0.0, 1.0]), [r1, r2], "with fixes"; fixes = [Fix("f", :typo, "x"; evidence = "e")])
    fg = Dict(figure_grades(Reproduction(id = "x", citation = "", doi = "", published = good, fixed = both)))
    @test fg["fig1"] == "A" && fg["fig2"] == "B"
end

@testset "post-comparison fixes are marked with ‡" begin
    ax = Axis("t", "-", :linear, 0, 1)
    r1 = ReferenceCurve("fig1_a", "f1", [0.5], [0.5], ax, ax, :curve, Dict{String, String}())
    bad = run_curves(_ -> ([0.0, 1.0], [1.0, 1.0]), [r1], "as published")
    fixed(fs) = run_curves(_ -> ([0.0, 1.0], [0.0, 1.0]), [r1], "with fixes"; fixes = fs)
    rep(fs) = Reproduction(id = "x", citation = "", doi = "", published = bad, fixed = fixed(fs))
    before = Fix("f", :typo, "x"; evidence = "e", basis = :paper, stage = :before_comparison)
    after = Fix("g", :units, "y"; evidence = "e", basis = :paper, stage = :after_comparison)
    target = Fix("h", :from_target, "z"; evidence = "e", basis = :figure, stage = :after_comparison)
    @test grade_label(rep([before])) == "B" && !post_comparison(rep([before]))
    @test grade_label(rep([before, after])) == "B‡" && post_comparison(rep([before, after]))
    @test grade_label(rep([target])) == "B†‡"
    @test Dict(figure_grades(rep([after])))["fig1"] == "B‡"
end

@testset "steep simulated segments are not clipped into a false match" begin
    ax = Axis("x", "-", :linear, 0, 1)
    ref = ReferenceCurve("c", "steep", [0.25], [0.5], ax, ax, :curve, Dict{String, String}())
    c = compare_curve(ref, [0.0, 1.0], [0.0, 10.0])        # line y = 10x passes far above (0.25, 0.5)
    @test isapprox(c.d[1], sqrt(0.25^2 + 0.5^2 - (0.25 + 5)^2 / 101); atol = 1e-12)   # true distance 0.199
    @test !passes(c)
    far = compare_curve(ref, [0.0, 1.0], [50.0, 50.0])     # whole curve far outside the box
    @test far.d[1] == 1.0
end

@testset "non-positive values on a log axis are gaps, not edge points" begin
    ax = Axis("x", "-", :linear, 0, 1)
    ref = ReferenceCurve("c", "log", [0.5], [sqrt(10.0)], ax, Axis("y", "-", :log, 1, 10), :curve, Dict{String, String}())
    @test compare_curve(ref, [0.0, 1.0], [-1.0, 100.0]).d == [1.0]          # no placeable segment: a complete miss
    c = compare_curve(ref, [0.0, 0.4, 1.0], [-1.0, 100.0, 100.0])
    @test !passes(c) && c.d[1] > 0.4
end


# ---- report-only diagnostics: coverage, reverse distance, band centre and edges -------------------
@testset "coverage of the reference and of the simulation" begin
    ax = Axis("t", "-", :linear, 0, 10)
    ay = Axis("y", "-", :linear, 0, 10)
    mk(xs; xaxis = ax) = ReferenceCurve("c", "f", Float64.(xs), Float64.(xs), xaxis, ay, :curve, Dict{String, String}())
    g = curve_diagnostics(mk([2.0, 5.0, 8.0]), 0.0:1.0:10.0, 0.0:1.0:10.0)
    @test g.coverage ≈ 0.6 && g.n_points == 3 && g.sim_coverage ≈ 1.0
    # a single point covers nothing; points beyond the axis are clamped to the whole axis
    @test curve_diagnostics(mk([5.0]), [0.0, 10.0], [0.0, 10.0]).coverage == 0.0
    @test curve_diagnostics(mk([-3.0, 14.0]), [0.0, 10.0], [0.0, 10.0]).coverage == 1.0
    # log x axis: span in decades (1 to 100 on 0.1 .. 1000 is 2 of 4 decades)
    lx = Axis("t", "-", :log, 0.1, 1000)
    @test curve_diagnostics(mk([1.0, 100.0]; xaxis = lx), [1.0, 10.0], [1.0, 10.0]).coverage ≈ 0.5
    # simulation defined on half of the axis only
    ref = mk([2.0, 8.0])
    @test curve_diagnostics(ref, [0.0, 5.0], [0.0, 5.0]).sim_coverage ≈ 0.5
    # a gap is not bridged: x = 4 and 5 are NaN, so [0, 3] and [6, 10] are defined
    xs = collect(0.0:1.0:10.0)
    ys = copy(xs); ys[5] = NaN; ys[6] = NaN
    @test curve_diagnostics(ref, xs, ys).sim_coverage ≈ 0.7
    # a lone finite point defines no interval; a simulation outside the axis is clipped to it
    @test curve_diagnostics(ref, [0.0, 5.0, 10.0], [NaN, 5.0, NaN]).sim_coverage == 0.0
    @test curve_diagnostics(ref, [-5.0, 20.0], [0.0, 10.0]).sim_coverage == 1.0
    # overlapping segments are counted once; unsorted x is fine
    @test curve_diagnostics(ref, [0.0, 6.0, 3.0, 9.0], [1.0, 1.0, 1.0, 1.0]).sim_coverage ≈ 0.9
    # non-positive values on a log y axis are gaps
    ly = Axis("y", "-", :log, 1, 10)
    rl = ReferenceCurve("c", "f", [2.0, 8.0], [2.0, 5.0], ax, ly, :curve, Dict{String, String}())
    @test curve_diagnostics(rl, [0.0, 5.0, 10.0], [3.0, 0.0, 3.0]).sim_coverage == 0.0
    @test curve_diagnostics(rl, [0.0, 5.0, 10.0], [-1.0, -1.0, -1.0]).sim_coverage == 0.0
    @test_throws DimensionMismatch curve_diagnostics(ref, [0.0, 1.0], [0.0])
end

@testset "reverse distance: the other side of a curve comparison" begin
    ax = Axis("t", "-", :linear, 0, 10)
    ay = Axis("y", "-", :linear, 0, 10)
    flat = ReferenceCurve("flat", "f", [0.0, 2.0, 4.0, 6.0, 8.0, 10.0], fill(5.0, 6), ax, ay, :curve, Dict{String, String}())
    # a simulation that crosses every reference point but wiggles between them passes the forward
    # check (every reference point lies on it) and is caught by the reverse distance
    xs = collect(0.0:0.01:10.0)
    wiggle = 5 .+ 3 .* sin.(pi .* xs ./ 2)
    c = compare_curve(flat, xs, wiggle)
    g = curve_diagnostics(flat, xs, wiggle)
    @test passes(c) && c.rms < 0.01
    @test g.rev_rms > 0.1 && g.rev_p95 > g.rev_rms && g.rev_n == length(xs)
    # the same simulation lying on the reference has reverse distance zero
    g0 = curve_diagnostics(flat, xs, fill(5.0, length(xs)))
    @test g0.rev_rms < 1e-12 && g0.rev_p95 < 1e-12
    # only simulated points inside the reference x range count
    ref = ReferenceCurve("r", "f", [2.0, 4.0, 6.0], [2.0, 4.0, 6.0], ax, ay, :curve, Dict{String, String}())
    gi = curve_diagnostics(ref, 0.0:1.0:10.0, 0.0:1.0:10.0)
    @test gi.rev_n == 5 && gi.rev_rms < 1e-12
    # a parallel line one unit above is 1/sqrt(2) of 10 away from the line y = x (interior points)
    gp = curve_diagnostics(ref, [3.0, 4.0, 5.0], [4.0, 5.0, 6.0])
    @test gp.rev_rms ≈ 1 / sqrt(2) / 10
    # capped at one axis box, like the forward distance
    @test curve_diagnostics(ref, [3.0, 4.0], [1e6, 1e6]).rev_rms == 1.0
    # gaps and unplaceable points are not used; nothing to measure gives NaN and rev_n = 0
    gg = curve_diagnostics(ref, [3.0, 4.0, 5.0], [3.0, NaN, 5.0])
    @test gg.rev_n == 2
    gn = curve_diagnostics(ref, [8.0, 9.0], [8.0, 9.0])
    @test gn.rev_n == 0 && isnan(gn.rev_rms) && isnan(gn.rev_p95)
    # a single reference point: distance to that point
    one = ReferenceCurve("o", "f", [5.0], [5.0], ax, ay, :curve, Dict{String, String}())
    @test curve_diagnostics(one, [5.0, 5.0], [6.0, 6.0]).rev_rms ≈ 0.1
    # the vertical metric is still a curve reference and gets a reverse distance
    vert = ReferenceCurve("v", "f", [2.0, 4.0, 6.0], [2.0, 4.0, 6.0], ax, ay, :curve, Dict("metric" => "vertical"))
    @test curve_diagnostics(vert, 0.0:1.0:10.0, 0.0:1.0:10.0).rev_n == 5
    # points and bars references have no reverse distance
    pts = ReferenceCurve("p", "f", [2.0, 4.0], [2.0, 4.0], ax, ay, :points, Dict{String, String}())
    gp = curve_diagnostics(pts, 0.0:1.0:10.0, 0.0:1.0:10.0)
    @test gp.rev_n == 0 && isnan(gp.rev_rms) && isnan(gp.band_centre_dev_rms)
end

@testset "band centre and simulated edges (report only)" begin
    ax = Axis("t", "-", :linear, 0, 10)
    ay = Axis("y", "-", :linear, 0, 10)
    band(xs, lo, hi; yaxis = ay) = ReferenceCurve("b", "f", Float64.(xs), (Float64.(lo) .+ Float64.(hi)) ./ 2, ax, yaxis, :band, Dict{String, String}(), Float64.(lo), Float64.(hi))
    ref = band([2.0, 4.0, 6.0], [1.0, 3.0, 5.0], [3.0, 5.0, 7.0])           # y in [x - 1, x + 1], half-width 1
    xs = collect(0.0:1.0:10.0)
    g = curve_diagnostics(ref, xs, xs)
    @test g.band_centre_dev_rms < 1e-12 && g.band_n == 3 && isnan(g.band_edge_rms) && isnan(g.band_edge_max)
    # inside the band but off the centre: the band score says 0, the centre deviation says 0.9
    off = xs .+ 0.9
    c = compare_curve(ref, xs, off)
    @test c.d == [0.0, 0.0, 0.0] && passes(c)
    @test curve_diagnostics(ref, xs, off).band_centre_dev_rms ≈ 0.9
    # outside by 2 half-widths on either side is 2, not capped
    @test curve_diagnostics(ref, xs, xs .- 2).band_centre_dev_rms ≈ 2.0
    # a gap leaves the reference x out instead of counting it
    gy = copy(xs); gy[4:6] .= NaN                                           # x = 3, 4, 5 undefined: x = 4 lies in the gap
    gg = curve_diagnostics(ref, xs, gy)
    @test gg.band_n == 2 && gg.band_centre_dev_rms < 1e-12
    # zero half-width reference points are left out; all zero gives NaN
    zero_hw = band([2.0, 4.0], [2.0, 3.0], [2.0, 5.0])
    @test curve_diagnostics(zero_hw, xs, xs .+ 0.5).band_centre_dev_rms ≈ 0.5      # only x = 4
    @test isnan(curve_diagnostics(band([2.0], [2.0], [2.0]), xs, xs).band_centre_dev_rms)
    # log axis: the centre is the geometric centre (displayed coordinate); 10 is the centre of [1, 100]
    ly = Axis("y", "-", :log, 0.01, 100)
    rl = band([5.0], [1.0], [100.0]; yaxis = ly)
    @test curve_diagnostics(rl, [0.0, 10.0], [10.0, 10.0]).band_centre_dev_rms < 1e-12
    @test curve_diagnostics(rl, [0.0, 10.0], [100.0, 100.0]).band_centre_dev_rms ≈ 1.0

    # simulated interval edges: compared with the reference edges at the reference x
    lo, hi = xs .- 1 .+ 0.2, xs .+ 1 .- 0.5                                  # lower edge 0.2 too high, upper 0.5 too low
    ge = curve_diagnostics(ref, xs, xs; sim_edges = (lo = lo, hi = hi))
    @test ge.band_edge_max ≈ 0.05 && ge.band_edge_rms ≈ sqrt((0.02^2 + 0.05^2) / 2)
    @test curve_diagnostics(ref, xs, xs; sim_edges = (lo = xs .- 1, hi = xs .+ 1)).band_edge_rms < 1e-12
    # edges are interpolated between simulated points, gaps are not bridged, errors are capped
    @test curve_diagnostics(ref, [0.0, 10.0], [0.0, 10.0]; sim_edges = (lo = [-1.0, 9.0], hi = [1.0, 11.0])).band_edge_rms < 1e-12
    @test isnan(curve_diagnostics(ref, xs, xs; sim_edges = (lo = fill(NaN, 11), hi = xs .+ 1)).band_edge_rms)
    @test curve_diagnostics(ref, xs, xs; sim_edges = (lo = fill(-1e6, 11), hi = fill(1e6, 11))).band_edge_max == 1.0
    @test_throws DimensionMismatch curve_diagnostics(ref, xs, xs; sim_edges = (lo = [1.0], hi = [2.0]))
    @test_throws ArgumentError curve_diagnostics(ref, xs, xs; sim_edges = (low = lo, high = hi))
    # edges are ignored for other kinds
    cur = ReferenceCurve("c", "f", [2.0], [2.0], ax, ay, :curve, Dict{String, String}())
    @test isnan(curve_diagnostics(cur, xs, xs; sim_edges = (lo = lo, hi = hi)).band_edge_rms)

    # the interpolation used for the diagnostics is the one of the band score
    refs = band([1.5, 3.2, 5.0, 7.7, 9.0], [0.0, 2.0, 4.0, 6.0, 7.0], [2.0, 4.0, 6.0, 8.0, 10.0])
    sx = [0.0, 1.0, 2.0, 4.0, 5.0, 8.0, 10.0]
    sy = [0.5, 1.5, NaN, 4.0, 5.5, 8.0, 12.0]
    s = QSPAtlas.sim_at_reference_x(refs, sx, sy)
    d = [isnan(si) ? 1.0 : min(max(lo_ - si, si - hi_, 0.0), 1.0) for (si, lo_, hi_) in zip(s, refs.ylo ./ 10, refs.yhi ./ 10)]
    @test compare_curve(refs, sx, sy).d ≈ d
end

@testset "fraction checked: the part of the simulated curve the reference compares" begin
    ax = Axis("t", "-", :linear, 0, 10)
    ay = Axis("y", "-", :linear, 0, 10)
    mk(xs; xaxis = ax, yaxis = ay, kind = :curve) = ReferenceCurve("c", "f", Float64.(xs), Float64.(xs), xaxis, yaxis, kind, Dict{String, String}())
    fc(ref, xs, ys) = curve_diagnostics(ref, xs, ys).fraction_checked
    # a straight line: arc length is proportional to x, so a reference on [2, 8] checks 60 %
    @test fc(mk([2.0, 5.0, 8.0]), [0.0, 10.0], [0.0, 10.0]) ≈ 0.6
    @test fc(mk([0.0, 10.0]), [0.0, 10.0], [0.0, 10.0]) ≈ 1.0
    # the part outside the axis box is not counted in either length: the line runs from -5 to 15
    @test fc(mk([2.0, 8.0]), [-5.0, 15.0], [-5.0, 15.0]) ≈ 0.6
    # a spike the reference does not reach: up the line x = 1 (0.9 long), then along y = 9 to x = 9
    # (0.8 long); a reference on [0, 5] checks 0.9 + 0.4 of 1.7
    @test fc(mk([0.0, 5.0]), [1.0, 1.0, 9.0], [0.0, 9.0, 9.0]) ≈ 13 / 17
    # a reference beyond the simulation, or a simulation drawn outside the box, leaves nothing to check
    @test fc(mk([9.0, 10.0]), [0.0, 5.0], [0.0, 5.0]) == 0.0
    @test isnan(fc(mk([2.0, 8.0]), [0.0, 10.0], [20.0, 30.0]))
    # a gap is not bridged: x = 4 and 5 are NaN, so [0, 3] and [6, 10] are drawn; [2, 8] checks 3 of 7
    xs = collect(0.0:1.0:10.0); ys = copy(xs); ys[5] = NaN; ys[6] = NaN
    @test fc(mk([2.0, 8.0]), xs, ys) ≈ 3 / 7
    # log axes are measured in decades: 1 to 100 on both axes, a reference on [10, 100] checks half
    lg = Axis("t", "-", :log, 1, 100)
    @test fc(mk([10.0, 100.0]; xaxis = lg, yaxis = lg), [1.0, 100.0], [1.0, 100.0]) ≈ 0.5
    # a band is checked over its own x range; points like a curve; bars have no meaning for it
    band = ReferenceCurve("b", "f", [2.0, 8.0], [5.0, 5.0], ax, ay, :band, Dict{String, String}(), [4.0, 4.0], [6.0, 6.0])
    @test fc(band, [0.0, 10.0], [0.0, 10.0]) ≈ 0.6
    @test fc(mk([2.0, 8.0]; kind = :points), [0.0, 10.0], [0.0, 10.0]) ≈ 0.6
    @test isnan(fc(mk([2.0, 8.0]; kind = :bars), [0.0, 10.0], [0.0, 10.0]))
end

@testset "reverse-distance criterion: opt-in, curve references only" begin
    ax = Axis("t", "-", :linear, 0, 10)
    ay = Axis("y", "-", :linear, 0, 10)
    # the reference ends at x = 2 and resumes at x = 8, y = 5 throughout; the simulation passes through
    # every reference point (forward error zero) but climbs to the top of the box between them
    rx = [0.0, 1.0, 2.0, 8.0, 9.0, 10.0]
    flat = ReferenceCurve("flat", "f", rx, fill(5.0, 6), ax, ay, :curve, Dict{String, String}())
    xs = collect(0.0:0.01:10.0)
    spike = [3 < x < 7 ? 10.0 : 5.0 for x in xs]
    clean = fill(5.0, length(xs))
    on = Tolerance(reverse_p95_max = 0.31)
    @test Tolerance().reverse_p95_max === nothing
    # off by default: the spike costs nothing
    c0 = compare_curve(flat, xs, spike)
    @test passes(c0) && !c0.failed_reverse && c0.rms < 1e-9
    @test curve_diagnostics(flat, xs, spike).rev_p95 > 0.31
    # on: the same curve fails, and says why; the forward numbers do not change
    c1 = compare_curve(flat, xs, spike; tolerance = on)
    @test !passes(c1) && c1.failed_reverse && c1.rms == c0.rms && c1.p95 == c0.p95 && c1.d == c0.d
    @test occursin("reverse", sprint(show, c1))
    # a simulation on the reference passes with the criterion on
    c2 = compare_curve(flat, xs, clean; tolerance = on)
    @test passes(c2) && !c2.failed_reverse
    # the maximum is the whole-trace P95 of curve_diagnostics: just above and just below
    p95 = curve_diagnostics(flat, xs, spike).rev_p95
    @test !passes(compare_curve(flat, xs, spike; tolerance = Tolerance(reverse_p95_max = p95 - 1e-6)))
    @test passes(compare_curve(flat, xs, spike; tolerance = Tolerance(reverse_p95_max = p95 + 1e-6)))
    # a curve that fails only the rms and p95 rule is not flagged as a reverse failure
    off_by = compare_curve(flat, xs, clean .+ 0.5; tolerance = on)
    @test !passes(off_by) && !off_by.failed_reverse
    # a curve that fails both is flagged for the reverse one too
    both = compare_curve(flat, xs, spike .+ 0.5; tolerance = on)
    @test !passes(both) && both.failed_reverse && both.rms > on.rms
    # the vertical metric is a curve reference and is subject to the criterion
    vert = ReferenceCurve("v", "f", rx, fill(5.0, 6), ax, ay, :curve, Dict("metric" => "vertical"))
    @test passes(compare_curve(vert, xs, spike)) && !passes(compare_curve(vert, xs, spike; tolerance = on))
    # bands, points and bars are never subject to it
    band = ReferenceCurve("b", "f", rx, fill(5.0, 6), ax, ay, :band, Dict{String, String}(), fill(4.0, 6), fill(6.0, 6))
    cb = compare_curve(band, xs, spike; tolerance = on)
    @test passes(cb) && !cb.failed_reverse && cb.d == compare_curve(band, xs, spike).d
    pts = ReferenceCurve("p", "f", rx, fill(5.0, 6), ax, ay, :points, Dict{String, String}())
    @test passes(compare_curve(pts, xs, spike; tolerance = on))
    bars = ReferenceCurve("r", "f", [1.0, 2.0], [5.0, 5.0], ax, ay, :bars, Dict{String, String}())
    @test passes(compare_curve(bars, [1.0, 2.0], [5.0, 5.0]; tolerance = on))
    # undefined reverse distance (no simulated point inside the reference x range) fails when on only
    mid = ReferenceCurve("m", "f", [2.0, 5.0, 8.0], fill(5.0, 3), ax, ay, :curve, Dict{String, String}())
    sparse = [0.0, 1.0, 9.0, 10.0]                         # the segment 1 to 9 spans the reference, no point inside
    @test passes(compare_curve(mid, sparse, fill(5.0, 4)))
    cu = compare_curve(mid, sparse, fill(5.0, 4); tolerance = on)
    @test !passes(cu) && cu.failed_reverse && cu.rms < 1e-9
    @test isnan(curve_diagnostics(mid, sparse, fill(5.0, 4)).rev_p95)
    narrow = ReferenceCurve("n", "f", [4.0, 6.0], [5.0, 5.0], ax, ay, :curve, Dict{String, String}())
    @test passes(compare_curve(narrow, [0.0, 10.0], [5.0, 5.0])) && !passes(compare_curve(narrow, [0.0, 10.0], [5.0, 5.0]; tolerance = on))
    # a simulation entirely below a log axis is a complete miss, and its reverse distance is undefined
    ly = Axis("y", "-", :log, 1, 10)
    rl = ReferenceCurve("l", "f", [2.0, 8.0], [2.0, 5.0], ax, ly, :curve, Dict{String, String}())
    @test !compare_curve(rl, [0.0, 10.0], [-1.0, -1.0]).failed_reverse
    @test compare_curve(rl, [0.0, 10.0], [-1.0, -1.0]; tolerance = on).failed_reverse

    # through run_curves and score.toml: grade, per-curve flag and the tolerance used
    refs = [flat, band]
    sim(_) = (xs, spike)
    plain = run_curves(sim, refs, "as published")
    strict = run_curves(sim, refs, "as published"; tolerance = on)
    @test QSPAtlas.matched(plain) && !QSPAtlas.matched(strict)
    @test grade(Reproduction(id = "x", citation = "", doi = "", published = plain)) == 'A'
    @test grade(Reproduction(id = "x", citation = "", doi = "", published = strict)) == 'C'
    # a per-curve tolerance replaces the whole Tolerance, so it decides on its own
    mixed = run_curves(sim, refs, "as published"; tolerance = on, tolerances = Dict("flat" => Tolerance()))
    @test QSPAtlas.matched(mixed)
    sp = read_score(write_score(tempname() * ".toml", Reproduction(id = "x", citation = "", doi = "", published = plain)))
    @test !any(haskey(c, "failed_reverse") || haskey(c, "tol_reverse_p95") for c in sp["published"]["curves"])
    @test all(haskey(c, "fraction_checked") for c in sp["published"]["curves"])
    ss = read_score(write_score(tempname() * ".toml", Reproduction(id = "x", citation = "", doi = "", published = strict)))
    cs = Dict(c["id"] => c for c in ss["published"]["curves"])
    @test cs["flat"]["failed_reverse"] === true && !cs["flat"]["passed"] && cs["flat"]["tol_reverse_p95"] == 0.31
    @test !haskey(cs["b"], "failed_reverse") && cs["b"]["passed"]
    @test ss["published"]["matched"] == false && haskey(cs["flat"], "rev_p95")
end

@testset "diagnostics ride along run_curves and write_score without changing a score" begin
    ax = Axis("t", "-", :linear, 0, 10)
    ay = Axis("y", "-", :linear, 0, 10)
    cur = ReferenceCurve("cur", "f", [2.0, 5.0, 8.0], [2.0, 5.0, 8.0], ax, ay, :curve, Dict{String, String}())
    bnd = ReferenceCurve("bnd", "f", [2.0, 4.0, 6.0], [2.0, 4.0, 6.0], ax, ay, :band, Dict{String, String}(), [1.0, 3.0, 5.0], [3.0, 5.0, 7.0])
    pts = ReferenceCurve("pts", "f", [2.0, 5.0], [2.0, 5.0], ax, ay, :points, Dict{String, String}())
    refs = [cur, bnd, pts]
    xs = collect(0.0:0.5:10.0)
    sim(ref) = (xs, xs .+ 0.3)
    plain = run_curves(sim, refs, "as published")
    with_edges = run_curves(refs, "as published") do ref
        ref.kind === :band ? (xs, xs .+ 0.3, (lo = xs .- 0.7, hi = xs .+ 1.3)) : (xs, xs .+ 0.3)
    end
    # the pass/fail and the distances do not depend on the extras
    for (a, b) in zip(plain.curves, with_edges.curves)
        @test a.comparison.passed == b.comparison.passed && a.comparison.d == b.comparison.d
        @test a.comparison.rms == compare_curve(a.ref, xs, xs .+ 0.3).rms
    end
    @test isnan(plain.curves[2].diagnostics.band_edge_rms)
    @test with_edges.curves[2].diagnostics.band_edge_rms ≈ 0.03       # both edges 0.3 off, in a 10-unit axis
    @test plain.curves[1].diagnostics.coverage ≈ 0.6 && plain.curves[1].diagnostics.rev_n == 13
    # a model that returns a 4-argument CurveResult constructor call still works, edges optional
    c = plain.curves[1]
    @test QSPAtlas.CurveResult(c.ref, c.x, c.y, c.comparison).diagnostics.coverage == c.diagnostics.coverage

    rep = Reproduction(id = "d", citation = "c", doi = "x", published = with_edges)
    path = write_score(tempname() * ".toml", rep)
    s = read_score(path)
    cs = Dict(c["id"] => c for c in s["published"]["curves"])
    @test Set(keys(cs)) == Set(["cur", "bnd", "pts"])
    for (c, r) in zip(with_edges.curves, refs)
        w = cs[r.id]
        @test w["passed"] == c.comparison.passed && w["rms"] == round(c.comparison.rms; sigdigits = 4)
        @test w["p95"] == round(c.comparison.p95; sigdigits = 4) && w["n_points"] == length(r.x)
        @test w["coverage"] == round(c.diagnostics.coverage; sigdigits = 4) && w["sim_coverage"] == 1.0
    end
    @test cs["cur"]["coverage"] ≈ 0.6 && haskey(cs["cur"], "rev_rms") && cs["cur"]["rev_n"] == 13
    @test !haskey(cs["cur"], "band_centre_dev_rms") && !haskey(cs["cur"], "band_edge_rms")
    @test haskey(cs["bnd"], "band_centre_dev_rms") && haskey(cs["bnd"], "band_edge_rms") && !haskey(cs["bnd"], "rev_rms")
    @test cs["bnd"]["band_centre_dev_rms"] ≈ 0.3 atol = 1e-3
    @test !haskey(cs["pts"], "rev_rms") && !haskey(cs["pts"], "band_centre_dev_rms")
    # without edges the edge fields are simply absent
    s2 = read_score(write_score(tempname() * ".toml", Reproduction(id = "d", citation = "c", doi = "x", published = plain)))
    @test !any(haskey(c, "band_edge_rms") for c in s2["published"]["curves"])
    # a simulation entirely below a log axis (a complete miss) still has diagnostics
    ly = Axis("y", "-", :log, 1, 10)
    rl = ReferenceCurve("l", "f", [2.0, 8.0], [2.0, 5.0], ax, ly, :curve, Dict{String, String}())
    miss = run_curves(_ -> ([0.0, 10.0], [-1.0, -1.0]), [rl], "as published")
    @test miss.curves[1].comparison.d == [1.0, 1.0] && miss.curves[1].diagnostics.sim_coverage == 0.0
    @test miss.curves[1].diagnostics.rev_n == 0
end
