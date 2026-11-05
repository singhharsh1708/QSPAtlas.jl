using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
import CairoMakie

Base.include(@__MODULE__, joinpath(@__DIR__, "model.jl"))  # Base.include: a bare Module() has no include
Base.include(@__MODULE__, joinpath(@__DIR__, "protocol.jl"))  # Base.include: a bare Module() has no include

const ID = "ouerdani_2015"
const CITATION = "Ouerdani A, et al. CPT Pharmacometrics Syst Pharmacol 2015;4:660-668"
const DOI = "10.1002/psp4.12001"
const PMCID = "PMC4716582"

# Figure 4 in its original layout: typical P (solid) and K (dashed), left mice, right patients.
function plot_fig4(path; fixes = FIXES)
    M = CairoMakie
    fig = M.Figure(size = (900, 420))
    panels = (("left", "Time (days)", (0, 75), (0, 250), "Predicted population tumor size P and capacity K"),
              ("right", "Time (months)", (0, 82), (40, 70), ""))
    for (j, (tag, xl, xlim, ylim, yl)) in enumerate(panels)
        ax = M.Axis(fig[1, j]; xlabel = xl, ylabel = yl, limits = (xlim, ylim))
        ref = ReferenceCurve("fig4_$(tag)_P", "", [0.0], [0.0], QSPAtlas.Axis("t", "", :linear, 0, 1), QSPAtlas.Axis("y", "", :linear, 0, 1), :curve, Dict{String, String}())
        x, P = simulate_curve(ref; fixes)
        _, Kc = simulate_curve(ReferenceCurve("fig4_$(tag)_K", "", [0.0], [0.0], ref.xaxis, ref.yaxis, :curve, Dict{String, String}()); fixes)
        M.lines!(ax, x, P; color = :black, linewidth = 2, label = "P")
        M.lines!(ax, x, Kc; color = :black, linewidth = 2, linestyle = :dash, label = "K")
        j == 2 && M.axislegend(ax; position = :rt)
    end
    M.save(path, fig)
    return path
end

function reproduce(; plots = true, figdir = joinpath(@__DIR__, "figures"))
    refs = sort(collect(values(load_references(joinpath(@__DIR__, "reference")))); by = r -> r.id)
    published = run_curves(ref -> simulate_curve(ref; fixes = PUBLISHED_FILLS), refs, "as published";
        fixes = PUBLISHED_FILLS, tolerance = TOLERANCE, tolerances = TOLERANCES)
    fixed = isempty(FIXES) ? nothing :
        run_curves(ref -> simulate_curve(ref; fixes = FIXES), refs, "with fixes";
            fixes = FIXES, tolerance = TOLERANCE, tolerances = TOLERANCES)
    rep = Reproduction(; id = ID, citation = CITATION, doi = DOI, pmcid = PMCID, published, fixed)
    write_score(joinpath(@__DIR__, "score.toml"), rep)
    if plots
        mkpath(figdir)
        plot_reproduction(rep, joinpath(figdir, "comparison.png"); ncols = 2)
        plot_fig4(joinpath(figdir, "fig4.png"))
    end
    return rep
end

(abspath(PROGRAM_FILE) == (@__FILE__)) && display(reproduce())
