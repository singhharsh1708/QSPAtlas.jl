# Draws the digitised points (green crosses) over the original figure images -> figures/digitise_check.png (Figure 3 on top, Figure 4 below)
#   julia --startup-file=no --project=tools/digitise models/fornari_2019/reference/check_overlay.jl
include(joinpath(@__DIR__, "digitise_run.jl"))
using FileIO

function overlay(imgkey, specs, cfg_dir, out; k = 4)
    img = loadimg(IMAGES[imgkey])
    H, W = size(img)
    big = repeat(img, inner = (k, k))
    for spec in specs
        spec["image"] == imgkey || continue
        ref = load_ref(joinpath(cfg_dir, spec["id"] * ".csv"))
        fxi = invlin(spec["xcal"]); fyi = invlog(spec["ycal"])
        for (x, y) in zip(ref...)
            c = fxi(x); r = fyi(y)
            for d in -3:3
                for (rr, cc) in ((r * k - k / 2 + d, c * k - k / 2), (r * k - k / 2, c * k - k / 2 + d))
                    ri, ci = round(Int, rr), round(Int, cc)
                    (1 <= ri <= H * k && 1 <= ci <= W * k) && (big[ri, ci] = RGB(0.0, 0.8, 0.0))
                end
            end
        end
    end
    save(out, big)
end
invlin(c) = v -> c[1] + (v - c[2]) * (c[3] - c[1]) / (c[4] - c[2])
invlog(c) = v -> c[1] + (log10(v) - log10(c[2])) * (c[3] - c[1]) / (log10(c[4]) - log10(c[2]))
function load_ref(path)
    xs = Float64[]; ys = Float64[]
    for l in eachline(path)
        (startswith(l, "#") || startswith(l, "x,y")) && continue
        a = split(l, ','); push!(xs, parse(Float64, a[1])); push!(ys, parse(Float64, a[2]))
    end
    return xs, ys
end
cfg = TOML.parsefile(joinpath(@__DIR__, "digitise.toml"))
mkpath(joinpath(@__DIR__, "..", "figures"))
tmp3 = tempname() * ".png"; tmp4 = tempname() * ".png"
overlay("g003", cfg["curve"], @__DIR__, tmp3)
overlay("g004", cfg["curve"], @__DIR__, tmp4)
a = load(tmp3); b = load(tmp4)
w = max(size(a, 2), size(b, 2))
pad(x) = size(x, 2) == w ? x : hcat(x, fill(eltype(x)(RGB(1, 1, 1)), size(x, 1), w - size(x, 2)))
save(joinpath(@__DIR__, "..", "figures", "digitise_check.png"), vcat(pad(a), pad(b)))
println("ok")
