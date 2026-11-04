# Regenerates reference/*.csv for lindauer_2017 Figure 4 from the pixel tips in digitise.toml.
#   julia --project=tools/digitise models/lindauer_2017/reference/digitise_bars.jl models/lindauer_2017/reference/digitise.toml
using TOML, Printf
cfg = TOML.parsefile(ARGS[1])
outdir = dirname(abspath(ARGS[1]))
src = cfg["source"]
calib(p1, v1, p2, v2) = p -> v1 + (p - p1) * (v2 - v1) / (p2 - p1)
rows = Tuple{String, Float64}[]   # (panel id, base value) for the base-case file
for P in cfg["panel"]
    f = calib(P["xcal"]...)
    for (mult, key, mtag) in (("5x", "tips_x5", "x5"), ("0.2x", "tips_x0p2", "x0p2"))
        id = "fig4_$(P["id"])_$(mtag)"
        open(joinpath(outdir, id * ".csv"), "w") do io
            println(io, "# id: ", id)
            println(io, "# figure: ", P["figure"], ", bar for parameter multiplied by ", mult)
            println(io, "# x: bar position [(index - 0.5)/2] linear 0 8")
            println(io, "# y: log10 tumour volume ratio, 26 wk / initial [-] linear ", P["axis_lo"], " ", P["axis_hi"])
            println(io, "# kind: bars")
            println(io, "# source: ", src)
            println(io, "# method: manual bar-tip picking, pixel columns in digitise.toml; xcal ", join(P["xcal"], " "), " (px1 v1 px2 v2); the bar values are tumour size changes at the end of the bar, not a curve over x")
            println(io, "# order: ", join(P["order"], ", "))
            println(io, "x,y")
            for (i, px) in enumerate(P[key])
                @printf(io, "%.4g,%.4f\n", (i - 0.5) / 2, f(px))
            end
        end
    end
end
# base case values from the Figure 4 caption (text, not digitised); x = 0.25 for 5 mg/kg and 0.75 for 0.5 mg/kg
open(joinpath(outdir, "fig4_base_case.csv"), "w") do io
    println(io, "# id: fig4_base_case")
    println(io, "# figure: Figure 4 caption, log ratio at the base parameter values: 5 mg/kg (x = 0.25) and 0.5 mg/kg (x = 0.75)")
    println(io, "# x: dose slot [(index - 0.5)/2] linear 0 8")
    println(io, "# y: log10 tumour volume ratio, 26 wk / initial [-] linear -3 0.5")
    println(io, "# kind: bars")
    println(io, "# source: Figure 4 caption text (PMC5270293): -1.29 for 5 mg/kg, -0.39 for 0.5 mg/kg")
    println(io, "# method: values read from the caption, exact to the printed digits (not digitised)")
    println(io, "x,y")
    for (x, P) in zip((0.25, 0.75), cfg["panel"])
        @printf(io, "%.4g,%.4f\n", x, P["base_value"])
    end
end
println("wrote Figure 4 reference CSVs to ", outdir)
