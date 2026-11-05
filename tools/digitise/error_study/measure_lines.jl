# Measure the drawn line width (perpendicular, in pixels) of the target curves in the real PMC
# images, so the synthetic renderer uses the same widths. Per column, coverage of a pixel is its
# colour distance from the white background divided by that of the line colour; summed down the
# column it gives the vertical line thickness, which is multiplied by cos(slope) (slope from the
# coverage centroids of neighbouring columns). Writes results/line_widths.csv.
#   julia --project=tools/digitise tools/digitise/error_study/measure_lines.jl
include(joinpath(@__DIR__, "..", "digitise.jl"))
using Printf

const W = RGB(1.0, 1.0, 1.0)
function width_profile(img, rows, cols, colour; maxdist = 0.35)
    cc = RGB(colour...)
    full = colourdist(cc, W)
    out = Float64[]; cen = Float64[]
    for c in cols
        s = 0.0; m = 0.0
        for r in rows
            p = img[r, c]
            # coverage only for pixels on the white-to-line-colour mixing line
            t = clamp(colourdist(p, W) / full, 0, 1)
            mix = RGB(1 - t * (1 - red(cc)), 1 - t * (1 - green(cc)), 1 - t * (1 - blue(cc)))
            colourdist(p, mix) > maxdist * t + 0.08 && continue
            s += t; m += t * r
        end
        push!(out, s); push!(cen, s > 0 ? m / s : NaN)
    end
    slope = [i in (1, length(cen)) ? NaN : (cen[i + 1] - cen[i - 1]) / 2 for i in eachindex(cen)]
    perp = [out[i] * cos(atan(slope[i])) for i in eachindex(out) if isfinite(slope[i])]
    return median(perp), median(out)
end

P = papers_root()
cases = [
  ("hutchinson_2016 Fig 2a blue", "hutchinson_2016/supplement/PSP4-5-636-g002.jpg", 60:118, 90:200, (0.0, 0.447, 0.741)),
  ("hutchinson_2016 Fig 2b blue", "hutchinson_2016/supplement/PSP4-5-636-g002.jpg", 60:118, 430:540, (0.0, 0.447, 0.741)),
  ("hutchinson_2016 Fig 6a black T", "hutchinson_2016/supplement/PSP4-5-636-g006.jpg", 30:85, 80:120, (0.0, 0.0, 0.0)),
  ("hutchinson_2016 Fig 6a magenta V", "hutchinson_2016/supplement/PSP4-5-636-g006.jpg", 30:85, 80:120, (1.0, 0.0, 1.0)),
  ("ouerdani_2015 Fig 4 right P solid", "ouerdani_2015/supplement/PSP4-4-660-g004.jpg", 140:200, 470:530, (0.0, 0.0, 0.0)),
  ("ouerdani_2015 Fig 4 left P solid", "ouerdani_2015/supplement/PSP4-4-660-g004.jpg", 190:235, 230:275, (0.0, 0.0, 0.0)),
  ("stroh_2019 Fig 3a orange", "stroh_2019/supplement/PSP4-8-676-g003.jpg", 100:200, 200:280, (0.85, 0.325, 0.098)),
  ("stroh_2019 Fig 3b purple", "stroh_2019/supplement/PSP4-8-676-g003.jpg", 140:248, 600:690, (0.494, 0.184, 0.556)),
  ("fornari_2019 Fig 3 CMP red", "fornari_2019/supplement/PSP4-8-858-g003.jpg", 40:120, 250:320, (0.85, 0.1, 0.15)),
]
open(joinpath(@__DIR__, "results", "line_widths.csv"), "w") do io
    println(io, "line,perpendicular_width_px,vertical_thickness_px")
    for (name, f, rows, cols, col) in cases
        img = loadimg(joinpath(P, f))
        wp, wv = width_profile(img, rows, cols, col)
        @printf(io, "%s,%.2f,%.2f\n", name, wp, wv)
        @printf("%-36s perp %.2f  vert %.2f\n", name, wp, wv)
    end
end
