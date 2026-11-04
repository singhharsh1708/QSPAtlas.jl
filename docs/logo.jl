# Draws the QSPAtlas logo: a concentration-time curve with digitised points sitting on it,
# on a dark rounded tile. Writes docs/src/assets/logo.svg (docs sidebar), logo.png (README)
# and favicon.ico (drawn at 64 px, converted with macOS `sips`).
#
#   julia --project=. docs/logo.jl
using CairoMakie

const NAVY = colorant"#16324f"
const TEAL = colorant"#2ec4b6"
const WHITE = colorant"#ffffff"

# rounded square from (0, 0) to (1, 1)
function rounded_square(r; n = 16)
    pts = Point2f[]
    for (cx, cy, a0) in ((1 - r, r, -π / 2), (1 - r, 1 - r, 0.0), (r, 1 - r, π / 2), (r, r, π))
        for a in range(a0, a0 + π / 2; length = n)
            push!(pts, Point2f(cx + r * cos(a), cy + r * sin(a)))
        end
    end
    return pts
end

function logo(px)
    fig = Figure(size = (px, px), backgroundcolor = :transparent, figure_padding = 0)
    ax = Axis(fig[1, 1]; backgroundcolor = :transparent, limits = (0, 1, 0, 1), aspect = DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    poly!(ax, rounded_square(0.2); color = NAVY)

    # axes of a small plot
    lw = px / 64
    lines!(ax, [0.2, 0.2, 0.82], [0.8, 0.2, 0.2]; color = (WHITE, 0.85), linewidth = lw, linecap = :round)

    # one-compartment oral absorption curve: rise, peak, decay
    t = range(0, 1; length = 200)
    c(x) = (exp(-2.2x) - exp(-11x)) / 0.62
    xs = 0.2 .+ 0.6 .* t
    ys = 0.2 .+ 0.52 .* c.(t)
    lines!(ax, xs, ys; color = TEAL, linewidth = 2.2lw, linecap = :round)

    # digitised points on the curve
    for s in (0.06, 0.18, 0.38, 0.62, 0.9)
        scatter!(ax, [0.2 + 0.6s], [0.2 + 0.52c(s)]; color = WHITE, markersize = px / 11, strokecolor = NAVY, strokewidth = lw / 1.5)
    end
    return fig
end

dir = joinpath(@__DIR__, "src", "assets")
mkpath(dir)
save(joinpath(dir, "logo.svg"), logo(256))
save(joinpath(dir, "logo.png"), logo(256))
save(joinpath(dir, "favicon.png"), logo(64); px_per_unit = 1)
if Sys.isapple()
    run(`sips -s format ico $(joinpath(dir, "favicon.png")) --out $(joinpath(dir, "favicon.ico"))`)
    rm(joinpath(dir, "favicon.png"))
end
println("wrote logo.svg, logo.png", Sys.isapple() ? ", favicon.ico" : ", favicon.png (convert it to favicon.ico)")
