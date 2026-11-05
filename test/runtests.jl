using QSPAtlas, Test

@testset "QSPAtlas" begin
    @testset "core" begin
        include("core.jl")
    end
    @testset "regimen" begin
        isfile(joinpath(@__DIR__, "regimen.jl")) && include("regimen.jl")
    end
    @testset "population" begin
        isfile(joinpath(@__DIR__, "population.jl")) && include("population.jl")
    end
    @testset "steadystate" begin
        isfile(joinpath(@__DIR__, "steadystate.jl")) && include("steadystate.jl")
    end
    @testset "summaries" begin
        isfile(joinpath(@__DIR__, "summaries.jl")) && include("summaries.jl")
    end

    # Every model's reproduce.jl runs in its own module. The grade and the pass/fail
    # of every curve must match the committed score.toml, so a change in a solver,
    # MTK, or the model code that moves a result fails CI instead of silently
    # changing the atlas.
    @testset "reproduce $id" for id in model_ids()
        dir = joinpath(models_dir(), id)
        scorefile = joinpath(dir, "score.toml")
        expected = isfile(scorefile) ? read_score(scorefile) : nothing
        m = Module(Symbol("Model_", id))
        Base.include(m, joinpath(dir, "reproduce.jl"))
        # reproduce() writes score.toml as a side effect; restore the committed file so
        # running the tests never changes tracked files
        committed = isfile(scorefile) ? read(scorefile) : nothing
        rep = try
            Base.invokelatest(m.reproduce; plots = false)
        finally
            committed === nothing || write(scorefile, committed)
        end
        @test rep isa Reproduction
        @test rep.id == id
        if expected === nothing
            @warn "$id has no committed score.toml; nothing to compare against"
        else
            @test string(grade(rep)) == expected["grade"]
            @test get(expected, "grade_label", expected["grade"]) == grade_label(rep)
            for (label, run) in (("published", rep.published), ("fixed", rep.fixed))
                exp_run = get(expected, label, Dict())
                # a run that existed must still exist, run, and score exactly the same curves:
                # dropping a failing target must fail CI
                @test (run === nothing) == isempty(exp_run)
                run === nothing && continue
                @test QSPAtlas.ran(run) == exp_run["ran"]
                want = Dict(c["id"] => c for c in get(exp_run, "curves", []))
                @test Set(c.ref.id for c in run.curves) == Set(keys(want))
                for c in run.curves
                    w = get(want, c.ref.id, nothing)
                    @test w !== nothing && w["passed"] == c.comparison.passed
                    # numbers must not drift either: score.toml stores 4 significant digits. Band
                    # targets are scored from Monte Carlo simulations. Julia 1.11 changed the tail
                    # branch of randn (log(rand) to log1p(-rand)), so about 7 % of normal draws differ
                    # from 1.10 and the simulated percentiles move by up to their Monte Carlo error
                    # (kept below 0.003 of the axis); they get that allowance. Pass/fail above must
                    # still match exactly, and the solver output itself is bitwise identical.
                    atol = c.ref.kind === :band ? 3e-3 : 1e-5
                    w === nothing || @test isapprox(c.comparison.rms, w["rms"]; rtol = 1e-3, atol)
                end
            end
        end
    end
end
