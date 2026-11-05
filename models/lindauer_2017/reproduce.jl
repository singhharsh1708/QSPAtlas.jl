using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

Base.include(@__MODULE__, joinpath(@__DIR__, "model.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "protocol.jl"))

const ID = "lindauer_2017"
const CITATION = "Lindauer A, et al. CPT Pharmacometrics Syst Pharmacol 2017;6:11-20"
const DOI = "10.1002/psp4.12130"
const PMCID = "PMC5270293"

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
        plot_reproduction(rep, joinpath(figdir, "comparison.png"))
        plot_figure4(Dict(x.id => x for x in refs), FIXES, joinpath(figdir, "figure4_regenerated.png"))
    end
    return rep
end

# Regenerated Figure 4 in the original layout: horizontal bars, white = 5x, black = 0.2x, axis crossing at the base value.
using CairoMakie: CairoMakie, Figure, barplot!, vlines!, Label, hidespines!
function plot_figure4(refs, fixes, path)
    s = setup_from_fixes(fixes)
    fig = Figure(size = (1000, 640))
    for (k, (dose, title, lo, hi)) in enumerate((("5mgkg", "5 mg/kg", -3.0, 0.5), ("0p5mgkg", "0.5 mg/kg", -2.0, 0.5)))
        order = String.(split(refs[string("fig4_", dose, "_x5")].meta["order"], ", "))
        base = cached_log_ratio(s, DOSE_MGKG[dose], "base", 1.0)
        ax = CairoMakie.Axis(fig[1, k]; title = "log ratio of tumour volume at 26 weeks vs initial ($title)", titlesize = 12,
            yticks = (1:16, order), yreversed = true, xaxisposition = :top, limits = (lo, hi, 0.3, 16.7), xlabel = "log10 ratio")
        for (j, name) in enumerate(order)
            for (mult, col, off) in ((5.0, :white, -0.2), (0.2, :black, 0.2))
                v = cached_log_ratio(s, DOSE_MGKG[dose], name, mult)
                isfinite(v) && CairoMakie.poly!(ax, CairoMakie.Rect(min(base, v), j + off - 0.2, abs(v - base), 0.4); color = col, strokecolor = :black, strokewidth = 1)
            end
        end
        vlines!(ax, [base]; color = :black)
        for (mult, tag, mk) in ((:x5, "x5", :circle), (:x0p2, "x0p2", :circle))
            r = refs["fig4_$(dose)_$(tag)"]
            CairoMakie.scatter!(ax, r.y, collect(1:16) .+ (tag == "x5" ? -0.2 : 0.2); color = tag == "x5" ? :red : :dodgerblue, markersize = 6)
        end
    end
    Label(fig[2, 1:2], "Regenerated (lindauer_2017, with fixes): white = 5x, black = 0.2x base value; red / blue dots = digitised bar tips from the paper";
        fontsize = 11)
    mkpath(dirname(path))
    CairoMakie.save(path, fig)
    return path
end

(abspath(PROGRAM_FILE) == (@__FILE__)) && display(reproduce())
