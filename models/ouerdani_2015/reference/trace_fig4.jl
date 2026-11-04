# Manual-rule digitisation of Figure 4 (PSP4-4-660-g004.jpg), both panels.
# Both P (solid) and K (dashed) are black, so colour tracing cannot separate them.
# Rule: per pixel column, take dark pixel runs (luminance < 0.6, rows 8..235 to skip
# frame/tick pixels), compute the darkness-weighted centroid row of each run, then
# assign runs to curves by position (P solid and continuous, K dashed):
#   LEFT : cols 51-103  P = lowest run on screen (largest row), K = runs above it
#          col 104, 106 skipped (curves cross / merged run)
#          col 105, 107-278  P = topmost run, K = all other runs
#   RIGHT: cols 311-350 single run = P (K is above the plot area)
#          cols 352-353 P = lowest run, K = runs above it (K enters the frame at col 351)
#          col 354 skipped (curves cross / merged run)
#          cols 355-537 P = topmost run, K = all other runs
# Calibration (pixels from tick marks, extra ticks checked in digitise_check.png):
#   left  x: col 48.0 = 0 d,  col 265.5 = 70 d    (check 50 d = col 203.5, 75 d = frame col 281)
#   left  y: row 239.5 = 0,   row 52.0 = 200 mm^3 (check 50/100/150 ticks, 250 = frame row 5)
#   right x: col 335.5 = 10 mo, col 506.0 = 70 mo (check 80 mo = col 534.5, 82 = frame col 540)
#   right y: row 43.5 = 65 mm, row 200.5 = 45 mm  (check 40/50/55/60/70 ticks)
# The image is not in the repository: it is read from $QSPATLAS_PAPERS/ouerdani_2015/supplement/ (tools/digitise/README.md).
# Run: julia --project=tools/digitise models/ouerdani_2015/reference/trace_fig4.jl [output folder]
# (the CSVs and digitise.toml are written to the output folder, default: this folder)
using Statistics, TOML, Printf
include(joinpath(@__DIR__, "..", "..", "..", "tools", "digitise", "digitise.jl"))

const OUTDIR = length(ARGS) > 0 ? abspath(ARGS[1]) : @__DIR__
mkpath(OUTDIR)
img = loadimg(paper_image("ouerdani_2015", "PSP4-4-660-g004.jpg"))

function column_runs(c; r1 = 8, r2 = 235, thr = 0.6)
    rs = [r for r in r1:r2 if lum(img[r, c]) < thr]
    return [(a, b, sum(r * (1 - lum(img[r, c])) for r in a:b) / sum(1 - lum(img[r, c]) for r in a:b)) for (a, b) in runs(rs)]
end

function pick(panel)
    P = Tuple{Int, Float64}[]; K = Tuple{Int, Float64}[]
    if panel == :left
        for c in 51:278
            rr = column_runs(c); isempty(rr) && continue
            cen = [x[3] for x in rr]
            if c <= 103
                push!(P, (c, cen[end])); for x in cen[1:end-1]; push!(K, (c, x)); end
            elseif c == 104
                push!(P, (c, cen[1]))   # merged run exactly at the crossing / peak of P
            elseif c == 106
                continue
            else
                push!(P, (c, cen[1])); for x in cen[2:end]; push!(K, (c, x)); end
            end
        end
    else
        for c in 311:537
            rr = column_runs(c); isempty(rr) && continue
            cen = [x[3] for x in rr]
            if c <= 350
                push!(P, (c, cen[1]))
            elseif c == 351
                push!(P, (c, cen[end]))
            elseif c <= 353
                push!(P, (c, cen[end])); for x in cen[1:end-1]; push!(K, (c, x)); end
            elseif c == 354
                continue
            else
                push!(P, (c, cen[1])); for x in cen[2:end]; push!(K, (c, x)); end
            end
        end
    end
    return P, K
end

cal = Dict(
    :left  => (fx = calib(48.0, 0.0, 265.5, 70.0, "linear"), fy = calib(239.5, 0.0, 52.0, 200.0, "linear")),
    :right => (fx = calib(335.5, 10.0, 506.0, 70.0, "linear"), fy = calib(43.5, 65.0, 200.5, 45.0, "linear")),
)
thin(v, n) = v[unique(round.(Int, range(1, length(v); length = min(n, length(v)))))]

function writecsv(id, figure, xaxis, yaxis, pts, panel; extra = "", n = 90)
    pts = thin(sort(pts), n)
    fx, fy = cal[panel].fx, cal[panel].fy
    open(joinpath(OUTDIR, id * ".csv"), "w") do io
        println(io, "# id: ", id)
        println(io, "# figure: ", figure)
        println(io, "# x: ", xaxis)
        println(io, "# y: ", yaxis)
        println(io, "# kind: curve")
        println(io, "# source: PSP4-4-660-g004.jpg (PMC4716582), Figure 4, authors' simulation with population parameters")
        println(io, "# method: manual rule-based pixel picking (black solid P and dashed K cannot be colour-traced); darkness-weighted run centroids per pixel column; rules and pixels in digitise.toml / trace_fig4.jl", extra)
        println(io, "x,y")
        for (c, r) in pts
            println(io, round(fx(c); sigdigits = 5), ",", round(fy(r); sigdigits = 5))
        end
    end
    return pts
end

out = Dict{String, Any}()
out["image"] = "ouerdani_2015/supplement/PSP4-4-660-g004.jpg"   # relative to $QSPATLAS_PAPERS
out["note"] = "Figure 4 digitised by manual rules, see trace_fig4.jl header for the assignment rules. Pixel coordinates (col,row) of every picked point are listed per curve before thinning."
out["calibration"] = Dict(
    "left" => Dict("xcal" => [48.0, 0.0, 265.5, 70.0], "ycal" => [239.5, 0.0, 52.0, 200.0], "xcheck" => "50 d at col 203.5, 75 d at frame col 281", "ycheck" => "ticks 50,100,150 at rows 192.4,145.6,98.5; 250 at frame row 5"),
    "right" => Dict("xcal" => [335.5, 10.0, 506.0, 70.0], "ycal" => [43.5, 65.0, 200.5, 45.0], "xcheck" => "80 mo at col 534.5, 82 mo at frame col 540", "ycheck" => "ticks 40,50,55,60,70 at rows 240,161.3,122,82.6,4.5"),
)
curves = Dict{String, Any}()
for (panel, ax, ay, tag, fig) in ((:left, "time [day] linear 0 75", "tumour volume [mm^3] linear 0 250", "left", "Figure 4 left, preclinical"),
                                  (:right, "time [month] linear 0 82", "tumour size SLD [mm] linear 40 70", "right", "Figure 4 right, clinical"))
    P, K = pick(panel)
    # P at t = 0: curve hugs the y axis (steep), first usable column is the one after the axis
    pp = writecsv("fig4_$(tag)_P", fig * ", tumour size P (solid)", ax, ay, P, panel)
    kp = writecsv("fig4_$(tag)_K", fig * ", carrying capacity K (dashed)", ax, ay, K, panel; n = 70)
    curves["fig4_$(tag)_P"] = Dict("npix_picked" => length(P), "n_csv" => length(pp), "cols" => first.(P), "rows" => round.(last.(P); digits = 2))
    curves["fig4_$(tag)_K"] = Dict("npix_picked" => length(K), "n_csv" => length(kp), "cols" => first.(K), "rows" => round.(last.(K); digits = 2))
end
out["curve"] = curves
open(joinpath(OUTDIR, "digitise.toml"), "w") do io
    TOML.print(io, out; sorted = true)
end
println("done")
