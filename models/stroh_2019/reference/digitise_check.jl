# (1) Independent tick check: converts tick pixels that were NOT used in the calibration into data values.
# (2) Overlay of all digitised points on the source images -> ../figures/digitise_check.png.
#   julia --project=tools/digitise models/stroh_2019/reference/digitise_check.jl
using TOML
include(joinpath(@__DIR__, "..", "..", "..", "tools", "digitise", "digitise.jl"))
cfg = TOML.parsefile(joinpath(@__DIR__, "digitise.toml"))
root = image_root(cfg, joinpath(@__DIR__, "digitise.toml"))
spec(id) = first(c for c in cfg["curve"] if c["id"] == id)
println("Tick check (pixel read from frame ticks not used for the calibration):")
# (curve id used for its calibration, axis, pixel, expected value)
for (id, ax, px, exp) in [
    ("fig3a_M1S1", :y, 116.5, 10.0), ("fig3a_M1S1", :y, 174.0, 1.0), ("fig3a_M1S1", :x, 144.5, 10.0), ("fig3a_M1S1", :x, 192.0, 15.0),
    ("fig3b_M1S2", :y, 126.0, 10.0), ("fig3b_M1S2", :y, 179.0, 1.0), ("fig3b_M1S2", :x, 512.5, 10.0),
    ("fig3c_M1S2", :y, 444.0, 10.0), ("fig3c_M1S2", :y, 395.5, 100.0), ("fig3c_M1S2", :x, 144.0, 10.0),
    ("fig4a_mAb", :y, 74.5, 10.0), ("fig4a_mAb", :y, 124.0, 0.1), ("fig4a_mAb", :x, 109.0, 20.0), ("fig4a_mAb", :x, 169.0, 40.0),
    ("fig4d_mAb", :y, 309.5, 1000.0), ("fig4d_mAb", :x, 111.0, 20.0), ("fig4e_mAb", :y, 284.3, 10.0), ("fig4e_mAb", :x, 348.0, 20.0)]
    s = spec(id)
    f = ax == :x ? calib(s["xcal"]..., s["xscale"]) : calib(s["ycal"]..., s["yscale"])
    v = f(px)
    println("  ", rpad(id, 12), ax, " px ", px, " -> ", round(v; sigdigits = 4), " (expected ", exp, ", rel. error ", round(abs(v - exp) / exp * 100; digits = 2), " %)")
end
function overlay(imgname, prefix; s = 2)
    img = loadimg(joinpath(root, imgname))
    big = repeat(img, inner = (s, s))
    for sp in cfg["curve"]
        startswith(sp["id"], prefix) || continue
        xc, yc = sp["xcal"], sp["ycal"]
        for l in eachline(joinpath(@__DIR__, sp["id"] * ".csv"))
            (startswith(l, "#") || startswith(l, "x,")) && continue
            x, y = parse.(Float64, split(l, ","))
            c = round(Int, (xc[1] + (x - xc[2]) * (xc[3] - xc[1]) / (xc[4] - xc[2])) * s)
            r = round(Int, (sp["yscale"] == "log" ? yc[1] + (log10(y) - log10(yc[2])) * (yc[3] - yc[1]) / (log10(yc[4]) - log10(yc[2])) :
                yc[1] + (y - yc[2]) * (yc[3] - yc[1]) / (yc[4] - yc[2])) * s)
            for dr in -2:2, dc in -2:2
                (abs(dr) == 2 || abs(dc) == 2) && 1 <= r + dr <= size(big, 1) && 1 <= c + dc <= size(big, 2) && (big[r + dr, c + dc] = RGB(0, 0, 0))
            end
        end
    end
    return big
end
a, b = overlay("PSP4-8-676-g003.jpg", "fig3"), overlay("PSP4-8-676-g004.jpg", "fig4")
W = max(size(a, 2), size(b, 2))
pad(m) = vcat(m, fill(RGB(1.0, 1.0, 1.0), 0, W)) |> x -> hcat(x, fill(RGB(1.0, 1.0, 1.0), size(x, 1), W - size(x, 2)))
save(joinpath(@__DIR__, "..", "figures", "digitise_check.png"), RGB.(vcat(pad(a), pad(b))))
println("wrote figures/digitise_check.png")
