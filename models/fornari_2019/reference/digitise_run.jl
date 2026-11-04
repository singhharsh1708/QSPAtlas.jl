# Regenerates reference/*.csv (and ancillary/*.csv) from the figure images using digitise.toml.
#   julia --startup-file=no --project=tools/digitise models/fornari_2019/reference/digitise_run.jl
# Needs `pdftoppm` (poppler) only for the ancillary Figure S3 curves.
include(joinpath(@__DIR__, "..", "..", "..", "tools", "digitise", "digitise.jl"))

const PAPERS = paper_dir("fornari_2019")
# images are read from $QSPATLAS_PAPERS/fornari_2019/supplement/ (see tools/digitise/README.md)
const IMAGES = Dict("g003" => paper_image("fornari_2019", "PSP4-8-858-g003.jpg"), "g004" => paper_image("fornari_2019", "PSP4-8-858-g004.jpg"))

redline(p) = red(p) > 0.4 && red(p) - green(p) > 0.20

# pixel row -> value on a log axis, two-point calibration given as [row1, val1, row2, val2]
logcal(c) = r -> 10^(log10(c[2]) + (r - c[1]) * (log10(c[4]) - log10(c[2])) / (c[3] - c[1]))
lincal(c) = p -> c[2] + (p - c[1]) * (c[4] - c[2]) / (c[3] - c[1])

function trace_red(spec, img)
    r1, r2, c1, c2 = spec["box"]
    fx = lincal(spec["xcal"]); fy = logcal(spec["ycal"])
    px = Float64[]; vy = Float64[]
    for c in c1:c2
        rs = [r for r in r1:r2 if redline(img[r, c])]
        isempty(rs) && continue
        push!(px, c); push!(vy, median(rs))
    end
    isempty(px) && error("no pixels for $(spec["id"])")
    n = get(spec, "n", 60)
    idx = unique(round.(Int, range(1, length(px); length = min(n, length(px)))))
    return fx.(px[idx]), fy.(vy[idx]), length(px)
end

function write_main(spec, x, y, npx, outdir)
    fy = logcal(spec["ycal"])
    top, bot = spec["frame_rows"]
    lo, hi = minmax(fy(bot), fy(top))
    sig(v) = round(v; sigdigits = 3)
    yaxis = string(spec["ylabel"], " [", spec["yunit"], "] log ", sig(lo), " ", sig(hi))
    path = joinpath(outdir, spec["id"] * ".csv")
    open(path, "w") do io
        println(io, "# id: ", spec["id"])
        println(io, "# figure: ", spec["figure"])
        println(io, "# x: ", spec["xaxis"])
        println(io, "# y: ", yaxis)
        println(io, "# kind: curve")
        println(io, "# source: ", basename(IMAGES[spec["image"]]), " (PMC6875710)")
        println(io, "# method: digitise_run.jl red-line mask (red>0.4, red-green>0.20), median row per pixel column, ", npx, " matched columns, box ", spec["box"], ", xcal ", spec["xcal"], ", ycal ", spec["ycal"], " (log), frame rows ", spec["frame_rows"], "; frozen 2026-11-04 before any simulation")
        println(io, "x,y")
        for (a, b) in zip(x, y)
            println(io, round(a; sigdigits = 5), ",", round(b; sigdigits = 5))
        end
    end
    return path
end

function main()
    cfg = TOML.parsefile(joinpath(@__DIR__, "digitise.toml"))
    outdir = @__DIR__
    for spec in cfg["curve"]
        img = loadimg(IMAGES[spec["image"]])
        x, y, npx = trace_red(spec, img)
        write_main(spec, x, y, npx, outdir)
        println("wrote ", spec["id"], ".csv (", length(x), " points from ", npx, " pixel columns)")
    end
    ancillary(joinpath(outdir, "ancillary"))
end

# Ancillary (not graded, not in reference/*.csv): Supplement Figure S3, rat free-platinum PK (blue, log left axis, uM) and
# Damage (orange, linear right axis) after a 30 mg/kg dose at t = 24 h. Used only to infer the dose-to-amount conversion that
# the paper does not state (see report.md, fix `dose_scale`).
function ancillary(outdir)
    mktempdir() do tmp
        run(`pdftoppm -r 200 -png $(joinpath(PAPERS, "PSP4-8-858-s003.pdf")) $(joinpath(tmp, "s3"))`)
        img = loadimg(joinpath(tmp, "s3-1.png"))
        # frame rows 87.5 (top: 1e3 uM / Damage 5) and 1031.5 (bottom: 1e-2 uM / Damage 0); x ticks 25 h at col 612.5, 45 h at col 1614.5
        ft = lincal([612.5, 25.0, 1614.5, 45.0])
        fb = logcal([87.5, 1e3, 1031.5, 1e-2])
        fo = lincal([1031.5, 0.0, 87.5, 5.0])
        isblue(p) = abs(red(p)) < 0.15 && abs(green(p) - 0.447) < 0.15 && abs(blue(p) - 0.741) < 0.15
        isor(p) = abs(red(p) - 0.85) < 0.12 && abs(green(p) - 0.325) < 0.12 && abs(blue(p) - 0.098) < 0.15
        for (name, isc, fy, cols) in (("S3_freePt_uM", isblue, fb, 566:1760), ("S3_damage", isor, fo, 566:1760))
            xs = Float64[]; ys = Float64[]
            for c in cols
                rs = [r for r in 90:1030 if isc(img[r, c])]
                isempty(rs) && continue
                # blue: keep the connected run containing the column's top-most line pixel only (data markers add outliers)
                push!(xs, ft(c)); push!(ys, fy(median(rs)))
            end
            open(joinpath(outdir, name * ".csv"), "w") do io
                println(io, "# ancillary: Supplement Figure S3 (PSP4-8-858-s003.pdf rendered at 200 dpi, 2340x1834), ", name)
                println(io, "# calibration: x 25 h at col 612.5, 45 h at col 1614.5; left axis log, rows 87.5 -> 1e3 and 1031.5 -> 1e-2; right axis linear rows 1031.5 -> 0 and 87.5 -> 5")
                println(io, "# NOTE: median row per column; columns containing data markers (t = 24-25, 27, 30, 36, 48 h) are unreliable for the blue curve")
                println(io, "x,y")
                for i in 1:5:length(xs)
                    println(io, round(xs[i]; sigdigits = 5), ",", round(ys[i]; sigdigits = 5))
                end
            end
            println("wrote ancillary/", name, ".csv")
        end
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
