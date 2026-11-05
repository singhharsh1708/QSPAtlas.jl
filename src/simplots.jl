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

# draw the quantile bands, lines and observed lines of a VPC into an axis
function draw_vpc!(ax, v::VPCResult; observed = nothing, logy = false)
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
    # the limit of quantification, where quantiles below it are left blank (censored VPC). Not on a
    # prediction-corrected axis: there every record has its own limit, LLOQ PRED_bin / PRED_ij, and
    # no single line is right, so none is drawn.
    if v.blq !== nothing && v.blq.handling === :censored && v.pc === :none && (!logy || v.blq.lloq > 0)
        CairoMakie.hlines!(ax, [v.blq.lloq]; color = :gray40, linestyle = :dot, linewidth = 1.5, label = "LLOQ")
    end
    axislegend(ax; position = :rt, framevisible = false, labelsize = 10)
    return ax
end

# the fraction of values below the limit: simulated band and median, observed points
function draw_blq_fraction!(ax, v::VPCResult)
    rows = sort(v.blq_fraction; by = r -> r.bin)
    x = [r.t_mid for r in rows]
    band!(ax, x, [r.lo for r in rows], [r.hi for r in rows]; color = (:steelblue, 0.25))
    lines!(ax, x, [r.median for r in rows]; color = :steelblue, linewidth = 2, label = "simulated")
    if v.observed !== nothing
        lines!(ax, x, [r.observed for r in rows]; color = :black, linewidth = 1.5)
        scatter!(ax, x, [r.observed for r in rows]; color = :black, markersize = 6, label = "observed")
    end
    axislegend(ax; position = :rt, framevisible = false, labelsize = 10)
    return ax
end

vpc_ylabel_suffix(v::VPCResult) = v.pc === :none ? "" : v.pc === :pc ? ", prediction-corrected" : ", prediction- and variability-corrected"

function vpc_title(v::VPCResult)
    return string(v.endpoint, v.stratum === nothing ? "" : string(" [", v.stratum, "]"), ": ", v.pc === :none ? "VPC" : v.pc === :pc ? "pcVPC" : "pvcVPC",
        ", ", v.replicates, " replicates, ", round(Int, 100 * v.pi_level), "% interval")
end

function vpc_panel!(pos, v::VPCResult; observed, logy, blq_fraction)
    grid = CairoMakie.GridLayout(pos)
    ax = MakieAxis(grid[1, 1]; title = vpc_title(v), xlabel = "time",
        ylabel = string(v.endpoint, " (", v.value, vpc_ylabel_suffix(v), ")"), yscale = logy ? log10 : identity)
    draw_vpc!(ax, v; observed, logy)
    if blq_fraction
        ax2 = MakieAxis(grid[2, 1]; xlabel = "time", ylabel = string("fraction below ", v.blq.lloq))
        draw_blq_fraction!(ax2, v)
        ylims!(ax2, 0, 1)
        CairoMakie.linkxaxes!(ax, ax2)
    end
    return ax
end

function check_blq_panel(v::VPCResult, blq_fraction)
    blq_fraction && v.blq_fraction === nothing && throw(ArgumentError("blq_fraction = true needs a VPC computed with blq = (lloq = ...)"))
    return nothing
end

"""
    plot_vpc(v::VPCResult; observed = nothing, logy = false, blq_fraction = false, path = nothing) -> Figure
    plot_vpc(v::StratifiedVPC; logy = false, blq_fraction = false, ncols = 3, path = nothing) -> Figure

Plot a [`vpc`](@ref) result. For every quantile a shaded band (the `pi_level` interval of the
replicate quantiles, `lo`..`hi`) and a line (their median); the median quantile is drawn in
orange, the others in blue. If the result has observed quantiles (`vpc(...; observed = data)`)
they are drawn as black lines with points. `observed` here may also be raw data, a vector of
`(id, time, dv)`, drawn as grey dots behind everything. Bins with an undefined quantile (NaN,
including quantiles below the limit in a censored VPC) leave gaps; with `blq.handling = :censored`
the limit is drawn as a dotted line, except for a prediction-corrected result: the correction
gives every record its own limit (`LLOQ * PRED_bin / PRED_ij`), so no single line would be
correct and none is drawn.

`blq_fraction = true` adds a panel below with the fraction of values below the limit of
quantification per bin (`VPCResult.blq_fraction`): the simulated median and interval, and the observed
fraction when the result has observed data; it needs a result computed with `blq`. A
prediction-corrected result is labelled as such.

For a [`StratifiedVPC`](@ref) every stratum gets its own panel (`ncols` per row, the stratum key in
the title; the raw `observed` points are not available there, the observed quantile lines are
drawn from the results).
"""
function plot_vpc(v::VPCResult; observed = nothing, logy::Bool = false, blq_fraction::Bool = false, path = nothing)
    check_blq_panel(v, blq_fraction)
    fig = Figure(size = (680, blq_fraction ? 700 : 460))
    vpc_panel!(fig[1, 1], v; observed, logy, blq_fraction)
    return finish_figure(fig, path)
end

function plot_vpc(sv::StratifiedVPC; logy::Bool = false, blq_fraction::Bool = false, ncols::Integer = 3, path = nothing)
    ncols >= 1 || throw(ArgumentError("ncols must be >= 1"))
    foreach(v -> check_blq_panel(v, blq_fraction), sv.results)
    n = length(sv)
    nc = min(n, ncols)
    nr = cld(n, nc)
    fig = Figure(size = (520 * nc, (blq_fraction ? 640 : 420) * nr))
    for (k, v) in enumerate(sv.results)
        vpc_panel!(fig[(k - 1) ÷ nc + 1, (k - 1) % nc + 1], v; observed = nothing, logy, blq_fraction)
    end
    return finish_figure(fig, path)
end
