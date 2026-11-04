using QSPAtlas, TOML, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "tools", "digitise", "papers.jl"))   # image_root: figure images from $QSPATLAS_PAPERS
# Check of the digitisation: calibration on ticks not used for calibration, plus overlay of the digitised points.
# Run from QSPAtlas.jl:  julia --project=. models/hutchinson_2016/reference/check_digitisation.jl
dir = joinpath(@__DIR__)
refs = load_references(dir)
cfg = TOML.parsefile(joinpath(dir, "digitise.toml"))
append!(cfg["curve"], TOML.parsefile(joinpath(dir, "digitise_fig6.toml"))["curve"])
root = image_root(cfg, joinpath(@__DIR__, "digitise.toml"))
# inverse calibration data -> pixel
inv1(v, (p1,v1,p2,v2)) = p1 + (v - v1)*(p2 - p1)/(v2 - v1)
fwd1(p, (p1,v1,p2,v2)) = v1 + (p - p1)*(v2 - v1)/(p2 - p1)
# extra ticks used only for checking: (curve id of the panel, axis, pixel, expected value)
checks = [
 ("fig2_a_T", :x, 192.0, 55), ("fig2_a_T", :x, 156.3, 50), ("fig2_a_T", :x, 227.5, 60), ("fig2_a_T", :y, 38.0, 300),
 ("fig2_b_V", :x, 539.0, 55), ("fig2_b_V", :y, 48.22, 600), ("fig2_b_V", :y, 19.46, 800),
 ("fig2_c_T_d04", :x, 192.0, 55), ("fig2_c_T_d04", :y, 224.0, 150),
 ("fig2_d_V_d04", :x, 539.0, 55), ("fig2_d_V_d04", :y, 224.0, 200), ("fig2_d_V_d04", :y, 205.0, 250), ("fig2_d_V_d04", :y, 243.22, 150),
 ("fig2_e_T_N2", :x, 192.0, 55), ("fig2_e_T_N2", :y, 374.0, 200), ("fig2_e_T_N2", :y, 397.0, 150),
 ("fig2_f_V_N2", :x, 539.0, 55), ("fig2_f_V_N2", :y, 385.31, 100),
 ("fig6_a_T", :x, 135.0, 50), ("fig6_a_T", :x, 183.0, 60), ("fig6_a_T", :x, 278.0, 80), ("fig6_a_T", :y, 6.83, 160), ("fig6_a_T", :y, 48.0, 80),
 ("fig6_b_T", :x, 442.6, 50), ("fig6_b_T", :x, 502.0, 60), ("fig6_b_T", :y, 6.82, 160), ("fig6_b_T", :y, 48.0, 80),
 ("fig6_c_T", :x, 129.0, 50), ("fig6_c_T", :x, 188.8, 60), ("fig6_c_T", :y, 133.13, 160), ("fig6_c_T", :y, 174.51, 80),
 ("fig6_d_T", :x, 442.6, 50), ("fig6_d_T", :x, 502.0, 60), ("fig6_d_T", :y, 133.0, 160), ("fig6_d_T", :y, 174.5, 80),
]
specs = Dict(s["id"] => s for s in cfg["curve"])
println("calibration check on ticks not used for calibration (data value from pixel vs printed label):")
maxerr = Dict(:x=>0.0, :y=>0.0)
for (id, ax, p, exp) in checks
    s = specs[id]
    v = fwd1(p, ax === :x ? s["xcal"] : s["ycal"])
    println(rpad(id, 14), ax, " pixel ", p, " -> ", round(v; digits=2), " expected ", exp)
end
function overlay!(fig, row, imgfile, prefix)
    img = CairoMakie.load(joinpath(root, imgfile))
    H, W = size(img)
    ax = CairoMakie.Axis(fig[row, 1]; yreversed = true, aspect = DataAspect(), title = imgfile * ": digitised points (red) on original")
    hidedecorations!(ax)
    image!(ax, (0.5, W + 0.5), (0.5, H + 0.5), permutedims(img))
    for (id, r) in refs
        startswith(id, prefix) || continue
        s = specs[id]
        px = inv1.(r.x, Ref(s["xcal"])); py = inv1.(r.y, Ref(s["ycal"]))
        scatter!(ax, px, py; color = :red, markersize = 4)
    end
    for (id, ax_, p, exp) in checks
        startswith(id, prefix) || continue
    end
end
fig = Figure(size = (1800, 2100))
overlay!(fig, 1, "PSP4-5-636-g002.jpg", "fig2")
overlay!(fig, 2, "PSP4-5-636-g006.jpg", "fig6")
rowsize!(fig.layout, 1, Relative(0.62))
save(joinpath(dir, "..", "figures", "digitise_check.png"), fig)
println("saved")
