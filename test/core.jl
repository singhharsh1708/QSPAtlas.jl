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

