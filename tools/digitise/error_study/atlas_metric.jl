# The atlas metric, loaded from the package source itself (src/reference.jl and src/compare.jl)
# without loading the simulation stack. `AtlasMetric.compare_curve` is QSPAtlas.compare_curve.
module AtlasMetric
using DelimitedFiles: readdlm, writedlm
using LinearAlgebra: norm
using Printf: @sprintf
using Statistics: mean, median, quantile
using TOML
const SRC = joinpath(@__DIR__, "..", "..", "..", "src")
include(joinpath(SRC, "reference.jl"))
include(joinpath(SRC, "compare.jl"))

"""Score digitised points (dx, dy) against the true curve (tx, ty) on the given axes."""
function score(dx, dy, tx, ty, xaxis::Axis, yaxis::Axis; id = "synthetic")
    ref = ReferenceCurve(id, id, Float64.(dx), Float64.(dy), xaxis, yaxis, :curve, Dict{String, String}())
    return compare_curve(ref, tx, ty)
end
end
