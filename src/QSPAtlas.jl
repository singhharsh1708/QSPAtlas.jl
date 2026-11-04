"""
    QSPAtlas

Reproduce, score and benchmark published PK/PD and QSP models.

Each model lives in `models/<author_year>/` and is reproduced by its `reproduce.jl`,
which compares simulations against points digitised from the paper's figures and
grades the result A-D (see [`grade`](@ref)).
"""
module QSPAtlas

using DelimitedFiles: readdlm, writedlm
using LinearAlgebra: norm
using Printf: @sprintf
using Statistics: mean, median, quantile
using TOML

include("reference.jl")
include("compare.jl")
include("score.jl")
include("dosing.jl")
include("regimen.jl")
include("population.jl")
include("simulate.jl")
include("summaries.jl")
include("plotting.jl")
include("simplots.jl")
include("benchmark.jl")

export Axis, ReferenceCurve, load_reference, load_references, write_reference
export Tolerance, CurveComparison, compare_curve, passes
export Fix, CurveResult, Run, Reproduction, grade, grade_label, from_target, figure_grades, run_curves, write_score, read_score
export models_dir, model_ids
export bolus_callback, repeated_doses
export Route, Dose, Regimen, DoseControl, DoseEvent, doses, resolve, dosing_callback
export ErrorModel, AdditiveError, ProportionalError, CombinedError, ExponentialError
export Subject, PopulationModel, SimResult, simulate, zero_re_prediction, virtual_population
export Reset, reset_times, SteadyStateOptions, SteadyStateError, steady_state_cycle, with_steady_state
export percentiles, vpc, VPCResult, spaghetti, bands, plot_vpc
export plot_reproduction
export benchmark_problem, BenchmarkSetup

"""
    models_dir()

Absolute path of the `models/` folder in this repository.
"""
models_dir() = normpath(joinpath(@__DIR__, "..", "models"))

"""
    model_ids()

Sorted ids (`author_year`) of all models that have a `reproduce.jl`.
"""
function model_ids()
    d = models_dir()
    return sort([m for m in readdir(d) if isfile(joinpath(d, m, "reproduce.jl"))])
end

end
