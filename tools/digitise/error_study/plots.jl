# Figures for the error study README (CairoMakie from the repository environment, read-only use):
#   julia --project=. tools/digitise/error_study/plots.jl <rms_tol> <p95_tol>
# Writes results/error_distributions.png and results/scores_vs_tolerance.png.
using CairoMakie, DelimitedFiles, Statistics, TOML

const RES = joinpath(@__DIR__, "results")
const MODELS = joinpath(@__DIR__, "..", "..", "..", "models")
readcsv(path) = (r = readdlm(path, ','; header = true); Dict(String(r[2][k]) => r[1][:, k] for k in eachindex(r[2])))
tr = parse(Float64, ARGS[1]); tp = parse(Float64, ARGS[2])

d = readcsv(joinpath(RES, "synthetic_runs.csv"))
real = readcsv(joinpath(RES, "real_figures.csv"))
panels = unique(d["panel"])
cols = Makie.wong_colors()
ecdf(v) = (s = sort(v); (s, (1:length(s)) ./ length(s)))

fig = Figure(size = (1100, 460))
for (j, (key, title, cur, new)) in enumerate((("rms", "RMS of point distances", 0.015, tr), ("p95", "95th percentile of point distances", 0.03, tp)))
    ax = Axis(fig[1, j]; xscale = log10, title, xlabel = "error (fraction of the axis box)", ylabel = "cumulative fraction of curves",
              xticks = ([1e-4, 1e-3, 1e-2, 1e-1], ["0.0001", "0.001", "0.01", "0.1"]))
    for (i, p) in enumerate(panels)
        k = (d["panel"] .== p) .& (d["traced"] .== true)
        x, y = ecdf(max.(Float64.(d[key][k]), 1e-4))
        lines!(ax, x, y; color = cols[mod1(i, 7)], linestyle = i > 7 ? :dash : :solid, label = p)
    end
    x, y = ecdf(Float64.(real[key]))
    lines!(ax, x, y; color = :black, linewidth = 3, label = "real: atlas vs thesis (34 curves)")
    vlines!(ax, [cur]; color = :grey40, linestyle = :dot, label = "current tolerance")
    vlines!(ax, [new]; color = :red, linestyle = :dash, label = "recommended tolerance")
    xlims!(ax, 1e-4, 0.3)
    j == 2 && Legend(fig[1, 3], ax; framevisible = false, labelsize = 11)
end
save(joinpath(RES, "error_distributions.png"), fig; px_per_unit = 1.5)

# committed per-curve scores (both runs) against the current and the recommended tolerance
fig2 = Figure(size = (720, 560))
ax = Axis(fig2[1, 1]; xscale = log10, yscale = log10, xlabel = "RMS (axis box)", ylabel = "P95 (axis box)",
          title = "Committed curve scores (models/*/score.toml)")
ms = sort(filter(x -> isfile(joinpath(MODELS, x, "score.toml")), readdir(MODELS)))
for (i, m) in enumerate(ms)
    s = TOML.parsefile(joinpath(MODELS, m, "score.toml"))
    for (run, mk) in (("published", :utriangle), ("fixed", :circle))
        cs = get(get(s, run, Dict()), "curves", [])
        isempty(cs) && continue
        scatter!(ax, max.([c["rms"] for c in cs], 1e-5), max.([c["p95"] for c in cs], 1e-5); color = cols[i], marker = mk, markersize = 8,
                 label = "$m ($run)")
    end
end
for (r, p, c, lab) in ((0.015, 0.03, :grey40, "current"), (tr, tp, :red, "recommended"))
    lines!(ax, [1e-5, r, r], [p, p, 1e-5]; color = c, linestyle = :dash, label = "$lab tolerance (pass = lower left)")
end
xlims!(ax, 3e-4, 1.2); ylims!(ax, 3e-4, 1.2)
Legend(fig2[1, 2], ax; framevisible = false, labelsize = 10)
save(joinpath(RES, "scores_vs_tolerance.png"), fig2; px_per_unit = 1.5)
println("saved")
