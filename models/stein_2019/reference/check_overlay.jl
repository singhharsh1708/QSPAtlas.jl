# Local check: draw the digitised band edges on the original figure (figures/digitise_check.png,
# never committed: it contains the paper's figure). Also reports the check against axis labels.
# Run: julia --startup-file=no --project=tools/digitise models/stein_2019/reference/check_overlay.jl
include(joinpath(@__DIR__, "trace_fig3a.jl"))
using FileIO, ImageCore
function overlay_refs(out; s = 3)
    big = repeat(img[1:270, :], inner = (s, s))
    colours = Dict("p90" => RGB(1, 0, 0), "p50" => RGB(0, 0.6, 0), "p10" => RGB(0.6, 0, 0.8))
    for (id, panel, band, _...) in TARGETS
        P = panels[panel]
        lines = split(read(joinpath(@__DIR__, id * ".csv"), String), '\n')
        for l in lines
            (isempty(l) || startswith(l, "#") || startswith(l, "x,")) && continue
            x, lo, hi = parse.(Float64, split(l, ','))
            col = P.xcal[1] + (x - P.xcal[2]) * (P.xcal[3] - P.xcal[1]) / (P.xcal[4] - P.xcal[2])
            for v in (lo, hi)
                row = P.ycal[1] + (log10(v) - log10(P.ycal[2])) * (P.ycal[3] - P.ycal[1]) / (log10(P.ycal[4]) - log10(P.ycal[2]))
                rr = round(Int, (row - 0.5) * s) + 1; cc = round(Int, (col - 0.5) * s) + 1
                (1 <= rr <= size(big, 1) - 2 && 1 <= cc <= size(big, 2) - 2) && (big[rr:(rr + 2), cc:(cc + 2)] .= colours[band])
            end
        end
    end
    save(out, RGB.(big))
end
mkpath(joinpath(@__DIR__, "..", "figures"))
overlay_refs(joinpath(@__DIR__, "..", "figures", "digitise_check.png"))
println("wrote figures/digitise_check.png")
