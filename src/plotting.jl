using CairoMakie: CairoMakie, Figure, Axis as MakieAxis, lines!, scatter!, axislegend, save, Label, xlims!, ylims!, band!

"""
    plot_reproduction(r::Reproduction, path; run = :best, ncols = 3)

Save a figure with one panel per reference curve: digitised reference points as
dots, the simulation as a line, and the pass/fail metrics in the panel title.
`run` picks `:published`, `:fixed`, or `:best` (fixed if present, else published).
"""
function plot_reproduction(r::Reproduction, path::AbstractString; run::Symbol = :best, ncols::Int = 3)
    R = run === :published ? r.published : run === :fixed ? r.fixed : something(r.fixed, r.published)
    R === nothing && throw(ArgumentError("$(r.id) has no run to plot"))
    n = max(length(R.curves), 1)
    nc = min(ncols, n)
    nr = cld(n, nc)
    fig = Figure(size = (380 * nc, 300 * nr + 50))
    Label(fig[0, 1:nc], string(r.id, ", ", R.label, ", grade ", grade_label(r)); fontsize = 16, font = :bold)
    for (k, c) in enumerate(R.curves)
        ref = c.ref
        ax = MakieAxis(fig[cld(k, nc), mod1(k, nc)];
            title = string(ref.id, "\n", c.comparison.passed ? "PASS" : "FAIL", @sprintf(" rms %.3f p95 %.3f", c.comparison.rms, c.comparison.p95)),
            titlesize = 11, xlabel = string(ref.xaxis.label, " [", ref.xaxis.unit, "]"),
            ylabel = string(ref.yaxis.label, " [", ref.yaxis.unit, "]"),
            yscale = ref.yaxis.scale === :log ? log10 : identity,
            xscale = ref.xaxis.scale === :log ? log10 : identity)
        # failed (non-finite) points stay as NaN so the line breaks there, as the scoring treats a gap
        ys = [isfinite(y) && (ref.yaxis.scale === :linear || y > 0) ? y : NaN for y in c.y]
        any(isfinite, ys) && lines!(ax, c.x, ys; color = :steelblue, linewidth = 2, label = "QSPAtlas")
        if ref.kind === :band
            # shaded interval; on a log axis only positive edges can be drawn
            ok = [i for i in eachindex(ref.x) if ref.yaxis.scale === :linear || ref.ylo[i] > 0]
            isempty(ok) || band!(ax, ref.x[ok], ref.ylo[ok], ref.yhi[ok]; color = (:black, 0.25), label = "digitised band")
        else
            scatter!(ax, ref.x, ref.y; color = :black, markersize = 7, label = "digitised")
        end
        xlims!(ax, ref.xaxis.lo, ref.xaxis.hi)
        ylims!(ax, ref.yaxis.lo, ref.yaxis.hi)
        k == 1 && axislegend(ax; position = :rt, labelsize = 9, framevisible = false)
    end
    mkpath(dirname(path))
    save(path, fig)
    return path
end
