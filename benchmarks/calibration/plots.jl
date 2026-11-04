# Figures and the estimates table for the hutchinson_2016 recalibration.
# Run from QSPAtlas.jl/:  julia --startup-file=no --project=benchmarks/calibration benchmarks/calibration/plots.jl
using StatsBase: Histogram, normalize
import StatsBase
using CairoMakie, CSV, DataFrames, TOML, Statistics, Distributions, LinearAlgebra, Printf
include(joinpath(@__DIR__, "common.jl"))
R = joinpath(@__DIR__, "results")
data = load_data(include_t0 = true)
pp = CSV.read(joinpath(R, "petab_predictions.csv"), DataFrame)
tp = CSV.read(joinpath(R, "turing_predictive.csv"), DataFrame)
mle = TOML.parsefile(joinpath(R, "petab_mle.toml"))
ts = CSV.read(joinpath(R, "turing_summary.csv"), DataFrame)
draws = CSV.read(joinpath(R, "turing_draws.csv"), DataFrame)
prof = CSV.read(joinpath(R, "petab_profiles.csv"), DataFrame)

cols = Dict("Control" => "#1f4e9c", "Bevacizumab" => "#2a9d4a", "Vanucizumab" => "#c8312b")
set_theme!(theme_light(); fontsize = 13)

# 1. fits
fig = Figure(size = (1100, 400))
for (j, g) in enumerate(GROUPS)
    ax = Axis(fig[1, j]; title = g, xlabel = "time (days after inoculation)", ylabel = j == 1 ? "mean tumour volume (mm³)" : "", limits = (36, 73, 0, 500))
    t = tp[tp.group .== g, :]
    band!(ax, t.time, t.q05, t.q95; color = (cols[g], 0.25), label = "Turing posterior 90% band (typical-individual curve)")
    lines!(ax, t.time, t.q50; color = cols[g], linewidth = 1.5, linestyle = :dot, label = "Turing posterior median")
    p = pp[pp.group .== g, :]
    lines!(ax, p.time, p.published; color = :black, linestyle = :dash, linewidth = 2, label = "published Table 1")
    lines!(ax, p.time, p.mle; color = cols[g], linewidth = 2.5, label = "PEtab MLE")
    d = data[data.group .== g, :]
    scatter!(ax, d.time, d.volume; color = :white, strokecolor = cols[g], strokewidth = 2, markersize = 11, label = "digitised group mean (Fig 1d)")
    j == 1 && axislegend(ax; position = :lt, framevisible = false, labelsize = 10)
end
save(joinpath(R, "fig_fits.png"), fig; px_per_unit = 2)

# 2. profiles
names = ["alpha_T", "alpha_V", "K", "Nmax", "delta_V_bev", "delta_V_van", "sigma"]
pubv = Dict("alpha_T" => PUB.alpha_T, "alpha_V" => PUB.alpha_V, "K" => PUB.K, "Nmax" => PUB.Nmax, "delta_V_bev" => PUB.delta_V_bev, "delta_V_van" => PUB.delta_V_van, "sigma" => 0.185)
nbest = mle["nllh"]
fig = Figure(size = (1200, 600))
for (i, nm) in enumerate(names)
    ax = Axis(fig[(i - 1) ÷ 4 + 1, (i - 1) % 4 + 1]; xlabel = nm, ylabel = "2 (nllh - nllh_min)", xscale = log10, limits = (nothing, (-0.5, 12)))
    d = sort(prof[prof.parameter .== "log10_" * nm, :], :value)
    lines!(ax, d.value, 2 .* (d.nllh .- nbest); color = :black, linewidth = 2)
    scatter!(ax, d.value, 2 .* (d.nllh .- nbest); color = :black, markersize = 5)
    hlines!(ax, [3.84]; color = :gray, linestyle = :dash)
    vlines!(ax, [pubv[nm]]; color = :red, linestyle = :dash)
    vlines!(ax, [mle[nm]]; color = :black, linestyle = :dot)
    r = ts[ts.parameter .== nm, :]
    nrow(r) == 1 && rangebars!(ax, [11.0], [r.q05[1]], [r.q95[1]]; direction = :x, color = :steelblue, linewidth = 4)
end
Label(fig[0, :], "PEtab profile likelihoods (black), 95% pointwise threshold 3.84 (grey), published value (red dashed), MLE (dotted), Turing 90% interval (blue bar)"; fontsize = 12)
save(joinpath(R, "fig_profiles.png"), fig; px_per_unit = 2)

# 3. posterior marginals vs prior
pm = (alpha_T = PUB.alpha_T, alpha_V = PUB.alpha_V, K = PUB.K, Nmax = PUB.Nmax, delta_V_bev = PUB.delta_V_bev, delta_V_van = PUB.delta_V_van, sigma = 0.2)
fig = Figure(size = (1200, 600))
for (i, nm) in enumerate(names)
    ax = Axis(fig[(i - 1) ÷ 4 + 1, (i - 1) % 4 + 1]; xlabel = nm, ylabel = "density (log10 axis)", xscale = log10)
    x = draws[!, Symbol(nm)]
    lx = log10.(x)
    h = StatsBase.fit(Histogram, lx, range(minimum(lx) - 0.05, maximum(lx) + 0.05; length = 40); closed = :left) |> h -> normalize(h; mode = :pdf)
    xs = 10 .^ (h.edges[1][1:(end - 1)] .+ step(h.edges[1]) / 2)
    stairs!(ax, xs, h.weights; color = :steelblue, linewidth = 2, step = :center)
    grid = 10 .^ range(log10(minimum(x)) - 0.7, log10(maximum(x)) + 0.7; length = 300)
    pd = nm == "Nmax" ? truncated(LogNormal(log(pm[Symbol(nm)]), 1.0); lower = 1.0) : LogNormal(log(pm[Symbol(nm)]), 1.0)
    lines!(ax, grid, pdf.(pd, grid) .* grid .* log(10); color = :gray, linestyle = :dash)   # prior density in log10 x
    vlines!(ax, [pubv[nm]]; color = :red, linestyle = :dash)
    vlines!(ax, [mle[nm]]; color = :black, linestyle = :dot)
end
Label(fig[0, :], "Turing posterior (blue) vs prior (grey dashed), published value (red dashed), PEtab MLE (dotted)"; fontsize = 12)
save(joinpath(R, "fig_posterior.png"), fig; px_per_unit = 2)

# 4. pairs: log-scale correlation heatmap and key scatter plots
C = cor(log.(Matrix(draws[:, Symbol.(names)])))
fig = Figure(size = (1200, 420))
ax = Axis(fig[1, 1]; title = "posterior correlation (log scale)", xticks = (1:7, names), yticks = (1:7, names), xticklabelrotation = pi / 4, yreversed = true, aspect = 1)
hm = heatmap!(ax, 1:7, 1:7, C; colormap = :balance, colorrange = (-1, 1))
for i in 1:7, j in 1:7
    text!(ax, j, i; text = @sprintf("%.2f", C[i, j]), align = (:center, :center), fontsize = 10)
end
Colorbar(fig[1, 2], hm)
for (k, (a, b)) in enumerate([("alpha_V", "delta_V_bev"), ("alpha_T", "alpha_V")])
    ax2 = Axis(fig[1, 2 + k]; xlabel = a, ylabel = b, xscale = log10, yscale = log10)
    scatter!(ax2, draws[!, Symbol(a)], draws[!, Symbol(b)]; markersize = 3, color = (:steelblue, 0.4))
    scatter!(ax2, [pubv[a]], [pubv[b]]; color = :red, marker = :star5, markersize = 16)
    scatter!(ax2, [mle[a]], [mle[b]]; color = :black, marker = :diamond, markersize = 12)
end
save(joinpath(R, "fig_pairs.png"), fig; px_per_unit = 2)

# 5. estimates table (toml)
est = Dict{String, Any}()
for nm in names
    r = ts[ts.parameter .== nm, :]
    est[nm] = Dict("published" => pubv[nm], "petab_mle" => mle[nm], "turing_median" => r.median[1], "turing_q05" => r.q05[1], "turing_q95" => r.q95[1],
        "rhat" => r.rhat[1], "ess_bulk" => r.ess_bulk[1])
end
open(joinpath(R, "estimates.toml"), "w") do io
    TOML.print(io, est)
end
println("plots done")
