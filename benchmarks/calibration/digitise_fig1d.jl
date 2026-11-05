# Digitise the MARKERS of Hutchinson 2016 Figure 1d (group mean tumour volume).
# Run: julia --startup-file=no --project=tools/digitise benchmarks/calibration/digitise_fig1d.jl   (from QSPAtlas.jl/)
# Method: classify pixels by colour (hue dominance), then for each measurement time take the vertical extent
# (bounding-box midpoint) of the coloured pixels in the +-2 px columns around the marker; the connecting line
# passes through the marker centre so it does not bias the midpoint. Overlay written for visual checking.
include(joinpath(@__DIR__, "..", "..", "tools", "digitise", "digitise.jl"))
using Statistics, FileIO

# The image is not in the repository: it is read from $QSPATLAS_PAPERS/hutchinson_2016/supplement/ (tools/digitise/README.md).
# Image provenance (same keys as the digitise.toml files): image_url = https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5192995/supplementaryFiles (member PSP4-5-636-g001.jpg),
# image_sha256 = 425e4178c358d7154e611c5cd0e27234adf46073fe7b85ebc13339486448a772
img = loadimg(paper_image("hutchinson_2016", "PSP4-5-636-g001.jpg"))

# calibration (pixels read with column_profile / row_profile, see notes.md in this folder)
# x: ticks 35 -> col 392, 75 -> col 662 ; y: 0 -> row 340, 800 -> row 199 (frame rows)
fx(c) = 35 + (c - 392) * 40 / (662 - 392)
fy(r) = (340 - r) * 800 / (340 - 199)
px(t) = 392 + (t - 35) * (662 - 392) / 40
py(v) = 340 - v * (340 - 199) / 800

# measurement times: twice a week from day 38 (10 per animal); read from the marker columns, see check below
times = [38, 42, 45, 49, 52, 56, 59, 63, 67, 71]

isblue(c) = blue(c) - red(c) > 0.3 && blue(c) - green(c) > 0.3
isblue_loose(c) = blue(c) - red(c) > 0.25 && blue(c) >= green(c) && !isgreen(c) && !isred(c)   # blue line/marker mixed with the green marker on top
isgreen(c) = green(c) - red(c) > 0.25 && green(c) - blue(c) > 0.25
isred(c) = red(c) - green(c) > 0.3 && red(c) - blue(c) > 0.3
groups = [("Control", isblue_loose), ("Bevacizumab", isgreen), ("Vanucizumab", isred)]

# legend box (excluded) and arrow rows (black, not matched anyway)
inlegend(r, c) = r <= 252 && c >= 545

# 1. check marker columns: for each group, columns where the coloured pixel count in the column is locally largest
for (g, f) in groups
    cnt = [count(r -> f(img[r, c]) && !inlegend(r, c), 200:339) for c in 393:661]
    cols = 393:661
    println(g, " column count maxima near expected marker columns:")
    for t in times
        c0 = round(Int, px(t))
        w = (c0 - 3):(c0 + 3)
        best = w[argmax([cnt[c - 392] for c in w])]
        println("  t=", t, " expected col ", c0, " best col ", best, " -> t=", round(fx(best); digits = 2))
    end
end

rows = NamedTuple[]
for (g, f) in groups
    for t in times
        c0 = round(Int, px(t)) - 1   # markers sit about 1 px left of the nominal tick position (check above)
        rs = Int[]
        for c in (c0 - 2):(c0 + 2), r in 200:339
            (f(img[r, c]) && !inlegend(r, c)) && push!(rs, r)
        end
        if isempty(rs)
            @warn "marker fully hidden, set to 70 mm^3 (Methods: mean size 70 mm^3 at day 38, groups randomised)" g t
            push!(rows, (; group = g, time = t, volume = 70.0, n = 0, span = 0))
            continue
        end
        # Control circles at 42 and 45 are mostly covered by the green markers/line, the colour rule picks up the
        # blue line leaving the marker: read by eye from a 12x zoom (circle top at row 323.2 and 321.5 + 2 px radius)
        if g == "Control" && t in (42, 45)
            v = t == 42 ? 85.0 : 95.0
            push!(rows, (; group = g, time = t, volume = v, n = 0, span = 0))
            continue
        end
        # keep the cluster closest to the median to drop crossing lines of the same colour
        m = median(rs)
        keep = filter(r -> abs(r - m) <= 6, rs)
        ymid = (minimum(keep) + maximum(keep)) / 2
        push!(rows, (; group = g, time = t, volume = round(fy(ymid); digits = 1), n = length(keep), span = maximum(keep) - minimum(keep)))
    end
end
for r in rows
    println(r)
end

open(joinpath(@__DIR__, "data", "hutchinson_fig1d.csv"), "w") do io
    println(io, "# Source: Hutchinson LG et al., CPT Pharmacometrics Syst Pharmacol 2016;5:636-645, Figure 1d (PSP4-5-636-g001.jpg, 667 px wide), mean tumour volume per group (n = 10 mice per group).")
    println(io, "# Method: marker centres digitised with tools/digitise/digitise.jl helpers via benchmarks/calibration/digitise_fig1d.jl (colour classification, bounding-box midpoint of marker in +-2 px columns).")
    println(io, "# Calibration: x ticks day 35 -> col 392, day 75 -> col 662; y frame row 340 = 0 mm^3, row 199 = 800 mm^3 (0.176 px per mm^3, 1 px = 5.7 mm^3). Measurement days 38,42,45,49,52,56,59,63,67,71 taken from the marker columns.")
    println(io, "# Digitisation error about +-1 px (about 6 mm^3). Hidden markers: day 38 Control and Bevacizumab are covered by the Vanucizumab triangle and set to 70 (Methods: mean size 70 mm^3 at day 38); Control days 42 and 45 are partly covered by the green markers and were read by eye from a 12x zoom (85 and 95, +-6).")
    println(io, "# Columns: group (Control = omalizumab isotype control, Bevacizumab, Vanucizumab), time [days after inoculation], volume [mm^3].")
    println(io, "group,time,volume")
    for r in rows
        println(io, r.group, ",", r.time, ",", r.volume)
    end
end

# overlay: 4x upscaled panel with digitised points as white-bordered crosses; also draws expected marker columns
crop = (150, 364, 340, 667)
r1, r2, c1, c2 = crop
sc = 4
H, W = (r2 - r1 + 1) * sc, (c2 - c1 + 1) * sc
big = [img[r1 + (i - 1) ÷ sc, c1 + (j - 1) ÷ sc] for i in 1:H, j in 1:W]
colfor = Dict("Control" => RGB(0.0, 0.0, 0.0), "Bevacizumab" => RGB(0.0, 0.0, 0.0), "Vanucizumab" => RGB(0.0, 0.0, 0.0))
function cross!(big, r, c, col; half = 9)   # hollow square around the digitised marker centre plus a centre dot
    big[r, c] = col
    for d in -half:half
        for (i, j) in ((r + d, c - half), (r + d, c + half), (r - half, c + d), (r + half, c + d))
            (1 <= i <= size(big, 1) && 1 <= j <= size(big, 2)) && (big[i, j] = col)
        end
    end
end
for r in rows
    c = (px(r.time) - c1) * sc + sc / 2 + 0.5
    rr = (py(r.volume) - r1) * sc + sc / 2 + 0.5
    cross!(big, round(Int, rr), round(Int, c), r.group == "Control" ? RGB(1.0, 0.0, 1.0) : r.group == "Bevacizumab" ? RGB(0.0, 0.0, 0.0) : RGB(0.0, 0.8, 0.8))
end
save(joinpath(@__DIR__, "results", "hutchinson_fig1d_overlay.png"), big)
println("done")
