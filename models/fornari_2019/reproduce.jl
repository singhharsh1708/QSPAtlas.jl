using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using CairoMakie: CairoMakie   # qualified use only: Makie's Axis would clash with QSPAtlas.Axis

Base.include(@__MODULE__, joinpath(@__DIR__, "model.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "protocol.jl"))

const ID = "fornari_2019"
const CITATION = "Fornari C, et al. CPT Pharmacometrics Syst Pharmacol 2019;8:858-868"
const DOI = "10.1002/psp4.12459"
const PMCID = "PMC6875710"

function reproduce(; plots = true, figdir = joinpath(@__DIR__, "figures"))
    refs = sort(collect(values(load_references(joinpath(@__DIR__, "reference")))); by = r -> r.id)
    published = run_curves(ref -> simulate_curve(ref; fixes = Fix[]), refs, "as published";
        tolerance = TOLERANCE, tolerances = TOLERANCES)
    fixed = isempty(FIXES) ? nothing :
        run_curves(ref -> simulate_curve(ref; fixes = FIXES), refs, "with fixes";
            fixes = FIXES, tolerance = TOLERANCE, tolerances = TOLERANCES)
    rep = Reproduction(; id = ID, citation = CITATION, doi = DOI, pmcid = PMCID, published, fixed)
    write_score(joinpath(@__DIR__, "score.toml"), rep)
    if plots
        plot_reproduction(rep, joinpath(figdir, "comparison.png"))
        plot_paper_figures(refs, figdir)
    end
    return rep
end

# Figures in the layout of the paper: Figure 3 (rat, 2 x 4) and Figure 4 a/b (human, 2 x 2).
# Black dots: digitised red model lines of the paper (graded). Red line: this reproduction (with fixes).
function plot_paper_figures(refs, figdir)
    M = CairoMakie
    byid = Dict(r.id => r for r in refs)
    fig = M.Figure(size = (1500, 650))
    layout = [("MPP", :MPP, "Multi-potent progenitors", "fig3_MPP"), ("CMP", :CMP, "Common myeloid progenitors", "fig3_CMP"),
        ("Neut", :Neut, "Neutrophils", "fig3_Neut"), ("Mono", :Mono, "Monocytes (not graded)", ""),
        ("MEP", :MEP, "Megakaryocyte-erythrocyte progenitors", "fig3_MEP"), ("Plt", :Plt, "Platelets", "fig3_Plt"),
        ("Ret", :Ret, "Reticulocytes (not graded)", ""), ("RBC", :RBC, "Red blood cells (not graded)", "")]
    order = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (2, 3), (2, 4)]
    for ((name, v, title, rid), (i, j)) in zip(layout, order)
        ax = M.Axis(fig[i, j]; title, xlabel = "Time (d)", ylabel = "1e9 cells / L", yscale = log10, titlesize = 13)
        x, y = series(:rat, v, FIXES)
        M.lines!(ax, x, y; color = :red, linewidth = 2, label = "this reproduction")
        if haskey(byid, rid)
            r = byid[rid]
            M.scatter!(ax, r.x, r.y; color = :black, markersize = 4, label = "paper (digitised)")
            M.ylims!(ax, r.yaxis.lo, r.yaxis.hi)
        end
        M.xlims!(ax, 0, 120)
        (i, j) == (1, 1) && M.axislegend(ax; position = :rt, framevisible = false, labelsize = 10)
    end
    M.Label(fig[0, 1:4], "Figure 3 regenerated: rat, carboplatin 40 mg/kg on day 1 of 8 repeated 14-day cycles"; fontsize = 16, font = :bold)
    M.save(joinpath(figdir, "fig3_regenerated.png"), fig)

    fig4 = M.Figure(size = (1000, 750))
    for (j, (case, panel, ttl)) in enumerate(((:human_a, "a", "(a) rat slope values"), (:human_b, "b", "(b) IC50-scaled slope values (Table 2)")))
        for (i, (v, lbl, thr, ylo, yhi)) in enumerate(((:Neut, "neut", (1.0, 0.5), 0.1, 100.0), (:Plt, "plt", (50.0, 25.0), 10.0, 1.0e4)))
            ax = M.Axis(fig4[i, j]; title = i == 1 ? ttl : "", xlabel = "Time (d)", ylabel = (v === :Neut ? "Neutrophils" : "Platelets") * " 1e9 cells / L", yscale = log10)
            x, y = series(case, v, FIXES)
            M.lines!(ax, x, y; color = :red, linewidth = 2, label = "this reproduction")
            r = byid["fig4$(panel)_$(lbl)"]
            M.scatter!(ax, r.x, r.y; color = :black, markersize = 4, label = "paper (digitised)")
            for th in thr
                M.hlines!(ax, th; color = :gray, linestyle = :dash)
            end
            M.xlims!(ax, 0, 42.4); M.ylims!(ax, ylo, yhi)
            (i, j) == (1, 1) && M.axislegend(ax; position = :rt, framevisible = false, labelsize = 10)
        end
    end
    M.Label(fig4[0, 1:2], "Figure 4 regenerated: human, single carboplatin dose (Zandvliet PK, AUC 5)"; fontsize = 16, font = :bold)
    M.save(joinpath(figdir, "fig4_regenerated.png"), fig4)
    return nothing
end

(abspath(PROGRAM_FILE) == (@__FILE__)) && display(reproduce())
