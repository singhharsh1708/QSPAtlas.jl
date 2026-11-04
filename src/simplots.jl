# Population plots (CairoMakie).
# QSPAtlas exports its own `Axis`, so the Makie one is always `MakieAxis` (imported in plotting.jl).

# Save the figure when a path is given; always return it.
function finish_figure(fig, path)
    if path !== nothing
        d = dirname(path)
        isempty(d) || mkpath(d)
        save(path, fig)
    end
    return fig
end

# Axis for one endpoint; on a log axis non-positive values are dropped by the callers (NaN).
function sim_axis(cell; endpoint, value, logy, title)
    return MakieAxis(cell; title, xlabel = "time", ylabel = string(endpoint, " (", value, ")"),
        yscale = logy ? log10 : identity)
end

positive_or_nan(v, logy) = logy && !(v > 0) ? NaN : v

"""
    spaghetti(r::SimResult, endpoint; value = :dv, replicate = 1, logy = false, path = nothing) -> Figure

One line per subject: `endpoint` against time. `replicate = 1` shows the first replicate of the
cohort; `replicate = nothing` overlays all of them. Non-finite values (failed solves) leave
gaps, they are not joined. `path` saves the figure (`png`, `pdf`, `svg` by extension).
"""
function spaghetti(r::SimResult, endpoint::Symbol; value::Symbol = :dv, replicate = 1, logy::Bool = false, path = nothing)
    check_value(value)
    items, rows = endpoint_items(r, endpoint, value)
    fig = Figure(size = (640, 420))
    ax = sim_axis(fig[1, 1]; endpoint, value, logy, title = string(endpoint, ": individual profiles"))
    byid = Dict{Tuple{Int, String}, Vector{Tuple{Float64, Float64}}}()
    for (o, (_, _, v)) in zip(rows, items)
        (replicate === nothing || o.replicate == replicate) || continue
        push!(get!(byid, (o.replicate, o.id), Tuple{Float64, Float64}[]), (o.time, positive_or_nan(v, logy)))
    end
    for key in sort!(collect(keys(byid)))
        pts = sort(byid[key]; by = first)
        lines!(ax, first.(pts), last.(pts); color = (:steelblue, 0.35), linewidth = 1)
    end
    return finish_figure(fig, path)
end

"""
    bands(r::SimResult, endpoint; q = (0.05, 0.5, 0.95), bins = :obstimes, value = :dv,
          logy = false, path = nothing) -> Figure

Percentile bands of a simulation: [`percentiles`](@ref) (pooled over replicates, same `bins`
and quantile definition) drawn as nested shaded bands between the lowest and the highest `q`,
the next lowest and next highest, and so on, and a line for an unpaired middle `q` (the median
for the default).
"""
function bands(r::SimResult, endpoint::Symbol; q = (0.05, 0.5, 0.95), bins = :obstimes, value::Symbol = :dv,
        logy::Bool = false, path = nothing, kwargs...)
    qs = sort(check_quantiles(q))
    rows = percentiles(r, endpoint; q = qs, bins, value, kwargs...)
    fig = Figure(size = (640, 420))
    ax = sim_axis(fig[1, 1]; endpoint, value, logy, title = string(endpoint, ": percentile bands"))
    draw_quantile_bands!(ax, rows, qs, r -> r.value, logy)
    return finish_figure(fig, path)
end

# per-q (x, y) vectors from rows
function qcurve(rows, qq, getter, logy)
    sel = sort([r for r in rows if r.q == qq]; by = r -> r.bin)
    return [r.t_mid for r in sel], [positive_or_nan(getter(r), logy) for r in sel]
end

# nested bands between q[i] and q[end+1-i]; an odd middle q as a line
function draw_quantile_bands!(ax, rows, qs, getter, logy)
    n = length(qs)
    for i in 1:(n ÷ 2)
        x, lo = qcurve(rows, qs[i], getter, logy)
        _, hi = qcurve(rows, qs[n + 1 - i], getter, logy)
        band!(ax, x, lo, hi; color = (:steelblue, 0.2 + 0.1 * (i - 1)), label = string(qs[i], " to ", qs[n + 1 - i]))
    end
    if isodd(n)
        x, y = qcurve(rows, qs[n ÷ 2 + 1], getter, logy)
        lines!(ax, x, y; color = :black, linewidth = 2, label = string("q = ", qs[n ÷ 2 + 1]))
    end
    axislegend(ax; position = :rt, framevisible = false, labelsize = 10)
    return ax
end

"""
    plot_vpc(v::VPCResult; observed = nothing, logy = false, path = nothing) -> Figure

Plot a [`vpc`](@ref) result. For every quantile a shaded band (the `pi_level` interval of the
replicate quantiles, `lo`..`hi`) and a line (their median); the median quantile is drawn in
orange, the others in blue. If the result has observed quantiles (`vpc(...; observed = data)`)
they are drawn as black lines with points. `observed` here may also be raw data, a vector of
`(id, time, dv)`, drawn as grey dots behind everything. Bins with an undefined quantile (NaN)
leave gaps.
"""
function plot_vpc(v::VPCResult; observed = nothing, logy::Bool = false, path = nothing)
    fig = Figure(size = (680, 460))
    ax = sim_axis(fig[1, 1]; endpoint = v.endpoint, value = v.value, logy,
        title = string(v.endpoint, ": VPC, ", v.replicates, " replicates, ", round(Int, 100 * v.pi_level), "% interval"))
    if observed !== nothing
        pts = [(Float64(o.time), positive_or_nan(Float64(o.dv), logy)) for o in observed]
        scatter!(ax, first.(pts), last.(pts); color = (:gray, 0.4), markersize = 5)
    end
    for qq in v.q
        col = qq == 0.5 ? :darkorange : :steelblue
        x, lo = qcurve(v.rows, qq, r -> r.lo, logy)
        _, hi = qcurve(v.rows, qq, r -> r.hi, logy)
        _, md = qcurve(v.rows, qq, r -> r.median, logy)
        band!(ax, x, lo, hi; color = (col, 0.25))
        lines!(ax, x, md; color = col, linewidth = 2, label = string("simulated q = ", qq))
        if v.observed !== nothing
            xo, yo = qcurve(v.observed, qq, r -> r.value, logy)
            lines!(ax, xo, yo; color = :black, linewidth = 1.5, linestyle = qq == 0.5 ? :solid : :dash)
            scatter!(ax, xo, yo; color = :black, markersize = 6)
        end
    end
    axislegend(ax; position = :rt, framevisible = false, labelsize = 10)
    return finish_figure(fig, path)
end
