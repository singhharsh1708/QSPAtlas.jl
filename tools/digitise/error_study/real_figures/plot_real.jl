# Overlay of the atlas reference points (from the PMC JPEGs) on the same curves taken from the theses.
# Only digitised data are drawn, never the thesis or paper images.
#   julia --project=. tools/digitise/error_study/real_figures/plot_real.jl     (root environment, for CairoMakie)
using CairoMakie, DelimitedFiles, Printf
include(joinpath(@__DIR__, "..", "atlas_metric.jl"))
include(joinpath(@__DIR__, "..", "..", "papers.jl"))
using .AtlasMetric: load_reference

const REPO = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
res = readdlm(joinpath(@__DIR__, "..", "results", "real_figures.csv"), ','; header = true)
score = Dict(String(r[2]) => (r[6], r[7]) for r in eachrow(res[1]))

panels = [
    ("ouerdani_2015", "Ouerdani Fig 4 left", ["fig4_left_P", "fig4_left_K"]),
    ("ouerdani_2015", "Ouerdani Fig 4 right", ["fig4_right_P", "fig4_right_K"]),
    ("hutchinson_2016", "Hutchinson Fig 2a", ["fig2_a_T"]),
    ("hutchinson_2016", "Hutchinson Fig 2b", ["fig2_b_V"]),
    [("hutchinson_2016", "Hutchinson Fig 2$p", ["fig2_$(p)_$(q)_$s" for s in ss]) for (p, q, ss) in (
        ("c", "T", ("d04", "d06", "d08", "d10", "d12")), ("d", "V", ("d04", "d06", "d08", "d10", "d12")),
        ("e", "T", ("N2", "N3", "N6", "N11", "N21")), ("f", "V", ("N2", "N3", "N6", "N11", "N21")))]...,
    [("hutchinson_2016", "Hutchinson Fig 6$p", ["fig6_$(p)_T", "fig6_$(p)_V"]) for p in ("a", "b", "c", "d")]...,
]

fig = Figure(size = (1500, 1050), fontsize = 12)
for (k, (model, title, ids)) in enumerate(panels)
    r, c = divrem(k - 1, 4)
    worst = maximum(score[id][1] for id in ids)
    ax = Axis(fig[r+1, c+1]; title = title, subtitle = @sprintf("max rms %.4f over %d curve(s)", worst, length(ids)),
        titlealign = :left, xgridvisible = false, ygridvisible = false)
    for id in ids
        ref = load_reference(joinpath(REPO, "models", model, "reference", id * ".csv"))
        t = readdlm(joinpath(papers_root(), model, "thesis", "curves", id * ".csv"), ',', Float64; skipstart = 1)
        lines!(ax, t[:, 1], t[:, 2]; color = (:gray25, 0.9), linewidth = 1.2, label = "thesis curve")
        scatter!(ax, ref.x, ref.y; color = (:darkorange2, 0.85), markersize = 4.5, label = "atlas reference (PMC JPEG)")
        xlims!(ax, ref.xaxis.lo, ref.xaxis.hi); ylims!(ax, ref.yaxis.lo, ref.yaxis.hi)
    end
    k == 1 && axislegend(ax; merge = true, position = :rt, framevisible = false, labelsize = 10)
end
Label(fig[0, :], "Atlas reference points against the same curves from the authors' theses (Ouerdani: 300 ppi raster; Hutchinson: vector)", fontsize = 14, halign = :left)
save(joinpath(@__DIR__, "..", "results", "real_figures.png"), fig; px_per_unit = 1.2)
println("wrote real_figures.png")
