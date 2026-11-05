# Score the atlas reference CSVs (digitised from the PMC JPEGs) against the same curves taken from the
# first authors' theses (higher-resolution raster or exact vector paths), with the atlas metric.
# The thesis curve plays the part of the "simulation": every atlas point gets its shortest distance to the
# thesis polyline in the atlas reference's own axis box (QSPAtlas.compare_curve, loaded from src/).
#
#   julia --project=tools/digitise tools/digitise/error_study/real_figures/compare_real.jl
#
# Thesis curves are read from $QSPATLAS_PAPERS/<id>/thesis/curves/<curve id>.csv (x,y in data units),
# written by ouerdani_thesis.jl and hutchinson_svg.py. Output: results/real_figures.csv.
using DelimitedFiles, Statistics, Printf
include(joinpath(@__DIR__, "..", "atlas_metric.jl"))
include(joinpath(@__DIR__, "..", "..", "papers.jl"))
using .AtlasMetric: load_reference, compare_curve, normalise

const REPO = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const RESULTS = joinpath(@__DIR__, "..", "results")
mkpath(RESULTS)

readxy(path) = (m = readdlm(path, ',', Float64; skipstart = 1); (m[:, 1], m[:, 2]))

# vertical offset of each reference point from the thesis curve, in axis-box units (positive: atlas above)
function signed_dy(ref, tx, ty)
    nx = normalise.(Ref(ref.xaxis), tx); ny = normalise.(Ref(ref.yaxis), ty)
    o = sortperm(nx); nx = nx[o]; ny = ny[o]
    out = Float64[]
    for (x, y) in zip(ref.x, ref.y)
        px = normalise(ref.xaxis, x); py = normalise(ref.yaxis, y)
        (px < nx[1] || px > nx[end]) && continue
        i = searchsortedlast(nx, px)
        i = clamp(i, 1, length(nx) - 1)
        s = nx[i+1] == nx[i] ? (ny[i] + ny[i+1]) / 2 : ny[i] + (px - nx[i]) / (nx[i+1] - nx[i]) * (ny[i+1] - ny[i])
        push!(out, py - s)
    end
    return out
end

sources = [
    ("ouerdani_2015", "thesis Fig 4 (HAL tel-01496174 PDF p96)", "raster 300 ppi", ["fig4_left_P", "fig4_left_K", "fig4_right_P", "fig4_right_K"]),
    ("hutchinson_2016", "thesis Fig 3.2 (ORA thesis PDF p85)", "vector",
        vcat("fig2_a_T", "fig2_b_V", ["fig2_$(p)_$(q)_$(s)" for (p, q, ss) in (("c", "T", ("d04", "d06", "d08", "d10", "d12")), ("d", "V", ("d04", "d06", "d08", "d10", "d12"))) for s in ss],
             ["fig2_$(p)_$(q)_$(s)" for (p, q) in (("e", "T"), ("f", "V")) for s in ("N2", "N3", "N6", "N11", "N21")])),
    ("hutchinson_2016", "thesis Fig 3.14 (ORA thesis PDF p111)", "vector", ["fig6_$(p)_$(q)" for p in ("a", "b", "c", "d") for q in ("T", "V")]),
]

rows = Any[]
pooled = Dict{String, Vector{Float64}}()
for (model, src, kind, ids) in sources, id in ids
    ref = load_reference(joinpath(REPO, "models", model, "reference", id * ".csv"))
    tx, ty = readxy(joinpath(papers_root(), model, "thesis", "curves", id * ".csv"))
    c = compare_curve(ref, tx, ty)
    dy = signed_dy(ref, tx, ty)
    push!(rows, (model, id, src, kind, length(ref.x), c.rms, c.p95, c.max, median(dy)))
    append!(get!(pooled, model * " " * split(id, '_')[1], Float64[]), c.d)
    @printf("%-16s %-14s n=%3d rms=%.4f p95=%.4f max=%.4f median_dy=%+.4f\n", model, id, length(ref.x), c.rms, c.p95, c.max, median(dy))
end
open(joinpath(RESULTS, "real_figures.csv"), "w") do io
    println(io, "model,curve_id,thesis_source,kind,n,rms,p95,max,median_signed_dy")
    for r in rows
        @printf(io, "%s,%s,%s,%s,%d,%.5f,%.5f,%.5f,%+.5f\n", r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9])
    end
end

# per figure: median and max of the per-curve rms and p95, and percentiles of all per-point distances
open(joinpath(RESULTS, "real_figures_summary.csv"), "w") do io
    println(io, "figure,curves,points,rms_median,rms_max,p95_median,p95_max,point_d_p50,point_d_p95,point_d_p99,point_d_max")
    for key in sort(collect(keys(pooled)))
        rs = [r for r in rows if r[1] * " " * split(r[2], '_')[1] == key]
        d = pooled[key]
        @printf(io, "%s,%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", key, length(rs), length(d), median(r[6] for r in rs), maximum(r[6] for r in rs),
            median(r[7] for r in rs), maximum(r[7] for r in rs), quantile(d, 0.5), quantile(d, 0.95), quantile(d, 0.99), maximum(d))
    end
end
