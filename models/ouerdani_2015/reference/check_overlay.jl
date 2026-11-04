# Overlay the digitised CSV points on the original figure (4x enlarged) and mark where
# the calibration maps extra ticks that were not used for calibration.
# Run: julia --project=tools/digitise models/ouerdani_2015/reference/check_overlay.jl
# The image is read from $QSPATLAS_PAPERS/ouerdani_2015/supplement/ (see tools/digitise/README.md).
using FileIO, ImageIO, Colors
include(joinpath(@__DIR__, "..", "..", "..", "tools", "digitise", "papers.jl"))
const ROOT = joinpath(@__DIR__, "..")
img = load(paper_image("ouerdani_2015", "PSP4-4-660-g004.jpg"))
k = 4
big = RGB{Float64}.(repeat(img, inner = (k, k)))
H, W = size(big)
function dot!(big, col, row, c; r = 2)
    i = round(Int, (row - 0.5) * k + 0.5); j = round(Int, (col - 0.5) * k + 0.5)
    for a in -r:r, b in -r:r
        (1 <= i + a <= size(big, 1) && 1 <= j + b <= size(big, 2)) && (big[i + a, j + b] = c)
    end
end
# data -> pixel (inverse of two-point calibration)
inv(p1, v1, p2, v2) = v -> p1 + (v - v1) * (p2 - p1) / (v2 - v1)
cx = Dict("left" => inv(48.0, 0.0, 265.5, 70.0), "right" => inv(335.5, 10.0, 506.0, 70.0))
cy = Dict("left" => inv(239.5, 0.0, 52.0, 200.0), "right" => inv(43.5, 65.0, 200.5, 45.0))
for (tag, cp, cK) in (("left", RGB(1, 0, 0), RGB(0, 0.6, 1)), ("right", RGB(1, 0, 0), RGB(0, 0.6, 1)))
    for (suffix, colr) in (("P", cp), ("K", cK))
        for line in eachline(joinpath(ROOT, "reference", "fig4_$(tag)_$(suffix).csv"))
            (startswith(line, "#") || startswith(line, "x")) && continue
            x, y = parse.(Float64, split(line, ','))
            dot!(big, cx[tag](x), cy[tag](y), colr; r = 1)
        end
    end
end
# extra ticks not used for calibration: green crosses at predicted position (x tick on axis row, y tick on axis col)
for (x, r) in ((50.0, 239.5), (75.0, 239.5)); dot!(big, cx["left"](x), r - 3, RGB(0, 0.8, 0); r = 2); end
for y in (50.0, 100.0, 150.0, 250.0); dot!(big, 52.0, cy["left"](y), RGB(0, 0.8, 0); r = 2); end
for (x, r) in ((80.0, 239.5), (82.0, 239.5)); dot!(big, cx["right"](x), r - 3, RGB(0, 0.8, 0); r = 2); end
for y in (40.0, 50.0, 55.0, 60.0, 70.0); dot!(big, 309.5, cy["right"](y), RGB(0, 0.8, 0); r = 2); end
save(joinpath(ROOT, "figures", "digitise_check.png"), big)
println("saved")
