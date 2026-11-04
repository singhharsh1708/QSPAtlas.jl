using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
import CairoMakie as CM

Base.include(@__MODULE__, joinpath(@__DIR__, "model.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "protocol.jl"))

const ID = "hutchinson_2016"
const CITATION = "Hutchinson LG, et al. CPT Pharmacometrics Syst Pharmacol 2016;5:636-645"
const DOI = "10.1002/psp4.12142"
const PMCID = "PMC5192995"

# MATLAB default line colours used in the paper's Figure 2
const COLS = [CM.RGBf(0, 0.447, 0.741), CM.RGBf(0.85, 0.325, 0.098), CM.RGBf(0.929, 0.694, 0.125), CM.RGBf(0.494, 0.184, 0.556), CM.RGBf(0.466, 0.674, 0.188)]

# Figure 2 in the paper's layout (a,c,e: T; b,d,f: V), with the digitised curves as dots
function plot_fig2(path, refs; fixes)
    fig = CM.Figure(size = (900, 750))
    layout = [("a", [("fig2_a_T", "")]), ("b", [("fig2_b_V", "")]),
              ("c", [("fig2_c_T_$k", "delta_V=" * string(DV[k])) for k in ("d04", "d06", "d08", "d10", "d12")]),
              ("d", [("fig2_d_V_$k", "delta_V=" * string(DV[k])) for k in ("d04", "d06", "d08", "d10", "d12")]),
              ("e", [("fig2_e_T_$k", "Nmax=" * string(Int(NM[k]))) for k in ("N2", "N3", "N6", "N11", "N21")]),
              ("f", [("fig2_f_V_$k", "Nmax=" * string(Int(NM[k]))) for k in ("N2", "N3", "N6", "N11", "N21")])]
    for (i, (panel, curves)) in enumerate(layout)
        r, c = fld1(i, 2), mod1(i, 2)
        isT = panel in ("a", "c", "e")
        ref1 = refs[curves[1][1]]
        ax = CM.Axis(fig[r, c]; xlabel = "time (days)", ylabel = isT ? "T (mm^3)" : "V (mm^3)", title = "($panel)",
            titlealign = :left, limits = ((ref1.xaxis.lo, ref1.xaxis.hi), (ref1.yaxis.lo, ref1.yaxis.hi)))
        for (j, (id, lab)) in enumerate(curves)
            x, y = simulate_curve(refs[id]; fixes)
            CM.lines!(ax, x, y; color = COLS[j], linewidth = 2, label = lab)
            CM.scatter!(ax, refs[id].x, refs[id].y; color = :black, markersize = 3)
        end
        panel in ("c", "e") && CM.axislegend(ax; position = :lt, labelsize = 10, framevisible = false)
    end
    CM.Label(fig[0, 1:2], "Hutchinson 2016 Figure 2 regenerated (lines: QSPAtlas with fixes; dots: digitised paper curves)"; fontsize = 13)
    mkpath(dirname(path)); CM.save(path, fig)
end

# Figure 6 in the paper's layout. Inputs as printed (T0 = 70 because Figure 6 gives none); the drug C is plotted
# as printed (Eq. 9), on its own right-hand axis, because its scale in the paper (0-5) is not what Eq. 9 gives.
function plot_fig6(path, refs; fixes)
    fig = CM.Figure(size = (1000, 600))
    for (i, panel) in enumerate(("a", "b", "c", "d"))
        p, _ = scenario("fig6_$(panel)_T")
        s = simulate_scenario(p; fixes, grid = collect(T_START:0.01:T_END))
        ax = CM.Axis(fig[fld1(i, 2), mod1(i, 2)]; xlabel = "time (days)", ylabel = "T and V (mm^3)", title = "($panel)", titlealign = :left,
            limits = ((35, 75), (0, 170)))
        CM.vspan!(ax, p.tnorm1, p.tnorm2; color = (:cyan, 0.3))
        p.tCoff > p.tCon && CM.vspan!(ax, p.tCon, p.tCoff; color = (:yellow, 0.3))
        CM.lines!(ax, s.t, s.T; color = :black, linewidth = 2)
        CM.lines!(ax, s.t, s.V; color = :magenta, linewidth = 2)
        for q in ("T", "V")
            r = get(refs, "fig6_$(panel)_$q", nothing)
            r === nothing || CM.scatter!(ax, r.x, r.y; color = q == "T" ? :gray40 : :purple, markersize = 3)
        end
        if p.tCoff > p.tCon
            ax2 = CM.Axis(fig[fld1(i, 2), mod1(i, 2)]; yaxisposition = :right, ylabel = "C (Eq. 9 as printed)", ylabelcolor = :blue,
                limits = ((35, 75), nothing), xticklabelsvisible = false, xticksvisible = false, xlabelvisible = false)
            CM.hidespines!(ax2, :t, :l, :b); CM.hidexdecorations!(ax2); CM.lines!(ax2, s.t, s.C; color = :blue, linestyle = :dash, linewidth = 1)
        end
    end
    CM.Label(fig[0, 1:2], "Hutchinson 2016 Figure 6 regenerated (black T, magenta V, dashed blue C; dots: digitised T and V)"; fontsize = 13)
    mkpath(dirname(path)); CM.save(path, fig)
end

"""
    reproduce(; plots = true)

Reproduce Figure 2 (22 curves) and Figure 6 (8 curves: T and V of panels a-d). All 30 were
frozen as targets together before any comparison, so all are graded (see report.md, Review).
"""
function reproduce(; plots = true, figdir = joinpath(@__DIR__, "figures"))
    refs = sort(collect(values(load_references(joinpath(@__DIR__, "reference")))); by = r -> r.id)
    allrefs = Dict(r.id => r for r in refs)
    published = run_curves(ref -> simulate_curve(ref; fixes = Fix[]), refs, "as published";
        tolerance = TOLERANCE, tolerances = TOLERANCES)
    fixed = isempty(FIXES) ? nothing :
        run_curves(ref -> simulate_curve(ref; fixes = FIXES), refs, "with fixes";
            fixes = FIXES, tolerance = TOLERANCE, tolerances = TOLERANCES)
    rep = Reproduction(; id = ID, citation = CITATION, doi = DOI, pmcid = PMCID, published, fixed)
    write_score(joinpath(@__DIR__, "score.toml"), rep)
    if plots
        plot_reproduction(rep, joinpath(figdir, "comparison.png"); ncols = 6)
        plot_fig2(joinpath(figdir, "fig2.png"), allrefs; fixes = FIXES)
        plot_fig6(joinpath(figdir, "fig6.png"), allrefs; fixes = FIXES)
    end
    return rep
end

(abspath(PROGRAM_FILE) == (@__FILE__)) && display(reproduce())
