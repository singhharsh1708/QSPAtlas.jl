using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
import CairoMakie

Base.include(@__MODULE__, joinpath(@__DIR__, "model.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "protocol.jl"))

const ID = "stroh_2019"
const CITATION = "Stroh M, et al. CPT Pharmacometrics Syst Pharmacol 2019;8:676-684"
const DOI = "10.1002/psp4.12448"
const PMCID = "PMC6765697"

const MATLAB = Dict("mAb00" => "#0072BD", "mAb" => "#0072BD", "M1S1" => "#D95319", "10x" => "#D95319", "M2S2" => "#77AC30", "10000x" => "#77AC30",
    "M1_0" => "#4DBEEE", "M2S1" => "#EDB120", "100x" => "#EDB120", "M1S2" => "#7E2F8E", "1000x" => "#7E2F8E")
const LABELS = Dict("mAb00" => "mAb0,0", "mAb" => "mAb0,0", "M1S1" => "Pb-Tx M1,S1", "M2S2" => "Pb-Tx M2,S2", "M1_0" => "Pb-Tx M1,0",
    "M2S1" => "Pb-Tx M2,S1", "M1S2" => "Pb-Tx M1,S2", "10x" => "10x Mask", "100x" => "100x Mask", "1000x" => "1000x Mask", "10000x" => "10000x Mask")

"Regenerate Figures 3 and 4 in the paper's layout: simulation lines (best run), digitised target as dots."
function paper_figures(rep, figdir)
    R = something(rep.fixed, rep.published)
    byid = Dict(c.ref.id => c for c in R.curves)
    mkpath(figdir)
    function panel!(ax, prefix, names)
        for n in names
            c = get(byid, prefix * n, nothing)
            c === nothing && continue
            ok = [(x, y) for (x, y) in zip(c.x, c.y) if isfinite(y) && (c.ref.yaxis.scale === :linear || y > 0)]
            CairoMakie.lines!(ax, first.(ok), last.(ok); color = MATLAB[n], linewidth = 2, label = LABELS[n])
            CairoMakie.scatter!(ax, c.ref.x, c.ref.y; color = :black, markersize = 2.5)
        end
    end
    f3 = CairoMakie.Figure(size = (900, 760))
    for (k, (p, dose, names)) in enumerate([("a", 3, ["mAb00", "M1S1", "M2S2", "M1_0"]), ("b", 5, ["mAb00", "M2S1", "M1S2"]), ("c", 10, ["mAb00", "M2S1", "M1S2"])])
        ax = CairoMakie.Axis(f3[cld(k, 2), mod1(k, 2)]; title = "($p) Dose = $dose mpk", xlabel = "Time (Days)", ylabel = "Total Drug (nM)", yscale = log10)
        panel!(ax, "fig3$(p)_", names)
        ref = byid["fig3$(p)_" * names[1]].ref
        CairoMakie.xlims!(ax, -1, 30); CairoMakie.ylims!(ax, ref.yaxis.lo, ref.yaxis.hi)
        CairoMakie.axislegend(ax; position = :rt, labelsize = 10)
    end
    CairoMakie.Label(f3[0, 1:2], "Stroh 2019 Figure 3 regenerated (lines: this model; dots: digitised published lines)"; fontsize = 13)
    CairoMakie.save(joinpath(figdir, "fig3_regenerated.png"), f3)
    names4 = ["mAb", "10x", "100x", "1000x", "10000x"]
    f4 = CairoMakie.Figure(size = (1200, 700))
    for (k, (p, ylab, log)) in enumerate([("a", "Total Drug (nM)", true), ("b", "Total Drug (nM)", true), ("c", "Total Drug (nM)", true), ("d", "Total Drug (nmol)", false), ("e", "Total Drug (nmol)", false)])
        ax = CairoMakie.Axis(f4[cld(k, 3), mod1(k, 3)]; title = "($p)", xlabel = "Time (Days)", ylabel = ylab, yscale = log ? log10 : identity)
        panel!(ax, "fig4$(p)_", names4)
        ref = byid["fig4$(p)_mAb"].ref
        CairoMakie.xlims!(ax, 0, 60); CairoMakie.ylims!(ax, ref.yaxis.lo, ref.yaxis.hi)
        k == 1 && CairoMakie.axislegend(ax; position = :rt, labelsize = 10)
    end
    CairoMakie.Label(f4[0, 1:3], "Stroh 2019 Figure 4 regenerated (human 4.5 mg/kg; lines: this model with fixes; dots: digitised published lines)"; fontsize = 13)
    CairoMakie.save(joinpath(figdir, "fig4_regenerated.png"), f4)
end

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
        paper_figures(rep, figdir)
    end
    return rep
end

(abspath(PROGRAM_FILE) == (@__FILE__)) && display(reproduce())
