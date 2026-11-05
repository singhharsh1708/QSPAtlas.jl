# Inventory of every digitised panel in the atlas: axis-box size in pixels (from the tick
# calibrations), axis scales, and the tracing rule used. Writes results/panels.csv.
#   julia --project=tools/digitise tools/digitise/error_study/panels.jl
using TOML, Printf
include(joinpath(@__DIR__, "atlas_metric.jl"))
using .AtlasMetric: load_reference, normalise

const MODELS = joinpath(@__DIR__, "..", "..", "..", "models")

# pixels spanned by the whole drawn axis (lo..hi) given a two-tick calibration [p1, v1, p2, v2]
axis_px(ax, cal) = abs(cal[3] - cal[1]) / abs(normalise(ax, cal[4]) - normalise(ax, cal[2]))

function calibrations()
    cal = Dict{String, Any}()
    for (m, files, rule) in (("hutchinson_2016", ["digitise.toml"], "colour_trace"),
                             ("hutchinson_2016", ["digitise_fig6.toml"], "dark_magenta_small"),
                             ("stroh_2019", ["digitise.toml"], "hue_track_log_markers"),
                             ("fornari_2019", ["digitise.toml"], "red_mask_log_bands"))
        for f in files, c in TOML.parsefile(joinpath(MODELS, m, "reference", f))["curve"]
            cal[c["id"]] = (c["xcal"], c["ycal"], rule, c["box"])
        end
    end
    o = TOML.parsefile(joinpath(MODELS, "ouerdani_2015", "reference", "digitise.toml"))["calibration"]
    for side in ("left", "right"), v in ("P", "K")
        cal["fig4_$(side)_$v"] = (o[side]["xcal"], o[side]["ycal"], "black_solid_dashed", nothing)
    end
    # Stein: constants of trace_fig3a.jl (panels dict)
    for p in ("p10", "p50", "p90")
        cal["fig3a_left_$p"] = ([82.5, 50.0, 275.5, 350.0], [90.5, 1e4, 192.5, 1e2], "band_edges", nothing)
        cal["fig3a_right_$p"] = ([406.0, 2.0, 634.0, 30.0], [89.0, 1e4, 185.0, 1e2], "band_edges", nothing)
    end
    return cal
end

function main()
    cal = calibrations()
    rows = String[]
    for m in sort(filter(d -> isdir(joinpath(MODELS, d, "reference")), readdir(MODELS))), f in sort(readdir(joinpath(MODELS, m, "reference")))
        endswith(f, ".csv") || continue
        r = load_reference(joinpath(MODELS, m, "reference", f))
        haskey(cal, r.id) || (println("skip ", m, " ", r.id, " (", r.kind, ")"); continue)
        xc, yc, rule, box = cal[r.id]
        push!(rows, join((m, r.id, rule, r.xaxis.scale, r.yaxis.scale, @sprintf("%.1f", axis_px(r.xaxis, xc)),
                          @sprintf("%.1f", axis_px(r.yaxis, yc)), length(r.x)), ","))
    end
    open(joinpath(@__DIR__, "results", "panels.csv"), "w") do io
        println(io, "model,curve_id,rule,xscale,yscale,x_axis_px,y_axis_px,n_points")
        foreach(l -> println(io, l), rows)
    end
    println(length(rows), " curves")
end
main()
