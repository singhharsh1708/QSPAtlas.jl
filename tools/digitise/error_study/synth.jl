# Synthetic digitising error study (Part A).
#
# For every panel type in the atlas, known curves are drawn at the panel's real pixel size with
# antialiasing (8x8 supersampling), with the line widths, colours, dash patterns, overlapping
# curves, shaded bands and markers of the real figures, saved as JPEG at the quality of the PMC
# files (IJG quality 92, 4:2:0 chroma subsampling for colour, greyscale for Ouerdani; see
# jpeg_quality.jl), reloaded, calibrated from tick marks found by sub-pixel centroids, traced with
# the same tracing code and settings the atlas uses for that panel type, and scored against the
# true curve with the atlas metric (QSPAtlas.compare_curve, loaded from src/).
#
#   julia --project=tools/digitise -t auto tools/digitise/error_study/synth.jl [reps] [seed]
#
# Output: results/synthetic_runs.csv (one row per traced curve), examples/ (a few synthetic JPEGs).
using Random, Statistics, Printf, FileIO, ImageIO, Colors, TOML

include(joinpath(@__DIR__, "atlas_metric.jl"))
using .AtlasMetric: Axis, score, normalise
# the atlas tracers, included unchanged: `trace_rects` (Hutchinson), `robust_trace` (Stroh), and
# through them `trace`, `calib`, `lum`, `chroma`, `colourdist`, `runs`, `loadimg` (tools/digitise/digitise.jl)
const MODELS = joinpath(@__DIR__, "..", "..", "..", "models")
include(joinpath(MODELS, "hutchinson_2016", "reference", "digitise_hutchinson.jl"))
include(joinpath(MODELS, "stroh_2019", "reference", "digitise_stroh.jl"))

const S = 8                      # supersampling factor per axis
const MARGIN = 10                # pixels around the axis box
const MATLAB = Dict("blue" => (0.0, 0.447, 0.741), "orange" => (0.85, 0.325, 0.098), "yellow" => (0.929, 0.694, 0.125),
                    "purple" => (0.494, 0.184, 0.556), "green" => (0.466, 0.674, 0.188), "lblue" => (0.301, 0.745, 0.933))
# colour targets the atlas uses for Hutchinson Figure 2 (digitise.toml), in MATLAB palette order
const HUT_TARGETS = Dict("blue" => [0.125, 0.375, 0.625], "orange" => [0.8, 0.35, 0.2], "yellow" => [0.9, 0.7, 0.25],
                         "purple" => [0.45, 0.2, 0.5], "green" => [0.45, 0.65, 0.25])

# ---------------------------------------------------------------- canvas and drawing
struct Canvas
    H::Int; W::Int                 # image size (pixels)
    XL::Float64; XR::Float64       # axis box edges, continuous pixel coordinates (pixel c spans c-0.5..c+0.5)
    YT::Float64; YB::Float64
    img::Matrix{RGB{Float64}}
end
function Canvas(wpx, hpx, rng)
    ox, oy = rand(rng), rand(rng)
    W = ceil(Int, wpx) + 2MARGIN; H = ceil(Int, hpx) + 2MARGIN
    XL = MARGIN + ox; YT = MARGIN + oy
    return Canvas(H, W, XL, XL + wpx, YT, YT + hpx, fill(RGB(1.0, 1.0, 1.0), H, W))
end
px_of(cv, u) = cv.XL + u * (cv.XR - cv.XL)          # normalised x -> pixel column coordinate
py_of(cv, v) = cv.YB - v * (cv.YB - cv.YT)          # normalised y -> pixel row coordinate
hires(cv) = falses(cv.H * S, cv.W * S)
hx(j) = 0.5 + (j - 0.5) / S                          # hi-res index -> continuous coordinate
hidx(X) = X * S - 0.5 * S + 0.5                       # inverse

function segd(px, py, ax, ay, bx, by)
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = L2 == 0 ? 0.0 : clamp(((px - ax) * dx + (py - ay) * dy) / L2, 0.0, 1.0)
    return hypot(px - ax - t * dx, py - ay - t * dy), t
end

"""Stroke a polyline (pixel coordinates) into mask `m`; optional dash pattern (on, off) in pixels; clipped to the axis box."""
function stroke!(m, cv, X, Y, w; dash = nothing, phase = 0.0, clip = true)
    r = w / 2
    s0 = 0.0
    for k in 1:(length(X) - 1)
        ax, ay, bx, by = X[k], Y[k], X[k + 1], Y[k + 1]
        L = hypot(bx - ax, by - ay)
        if isfinite(L)
            j1 = max(1, floor(Int, hidx(min(ax, bx) - r))); j2 = min(size(m, 2), ceil(Int, hidx(max(ax, bx) + r)))
            i1 = max(1, floor(Int, hidx(min(ay, by) - r))); i2 = min(size(m, 1), ceil(Int, hidx(max(ay, by) + r)))
            for j in j1:j2, i in i1:i2
                Xp, Yp = hx(j), hx(i)
                clip && !(cv.XL <= Xp <= cv.XR && cv.YT <= Yp <= cv.YB) && continue
                d, t = segd(Xp, Yp, ax, ay, bx, by)
                d > r && continue
                if dash !== nothing
                    s = s0 + t * L + phase
                    mod(s, dash[1] + dash[2]) > dash[1] && continue
                end
                m[i, j] = true
            end
            s0 += L
        end
    end
    return m
end

function ring!(m, cv, cx, cy, rad, w; shape = :circle)
    j1 = max(1, floor(Int, hidx(cx - rad - w))); j2 = min(size(m, 2), ceil(Int, hidx(cx + rad + w)))
    i1 = max(1, floor(Int, hidx(cy - rad - w))); i2 = min(size(m, 1), ceil(Int, hidx(cy + rad + w)))
    for j in j1:j2, i in i1:i2
        (cv.XL <= hx(j) <= cv.XR && cv.YT <= hx(i) <= cv.YB) || continue      # clipped to the axis box
        dx, dy = hx(j) - cx, hx(i) - cy
        d = shape === :circle ? hypot(dx, dy) : shape === :diamond ? (abs(dx) + abs(dy)) / sqrt(2) : max(abs(dx), abs(dy))
        abs(d - rad) <= w / 2 && (m[i, j] = true)
    end
    return m
end

"""Fill between two pixel-row functions lo(X) >= hi(X) (rows grow downwards) over columns X1..X2, clipped to the axis box."""
function fill_between!(m, cv, X1, X2, top, bot)
    for j in 1:size(m, 2)
        Xp = hx(j)
        (X1 <= Xp <= X2 && cv.XL <= Xp <= cv.XR) || continue
        t, b = top(Xp), bot(Xp)
        for i in 1:size(m, 1)
            Yp = hx(i)
            (t <= Yp <= b && cv.YT <= Yp <= cv.YB) && (m[i, j] = true)
        end
    end
    return m
end

"""Blend a hi-res mask onto the image with colour c and opacity a (coverage = fraction of sub-pixels)."""
function composite!(cv, m, c; a = 1.0)
    col = RGB(c...)
    for j in 1:cv.W, i in 1:cv.H
        n = 0
        @inbounds for jj in ((j - 1) * S + 1):(j * S), ii in ((i - 1) * S + 1):(i * S)
            n += m[ii, jj]
        end
        n == 0 && continue
        α = a * n / S^2
        p = cv.img[i, j]
        cv.img[i, j] = RGB((1 - α) * red(p) + α * red(col), (1 - α) * green(p) + α * green(col), (1 - α) * blue(p) + α * blue(col))
    end
end

"""Frame (1 px) and outward tick marks (below the x axis, left of the y axis) at the tick positions.
Outward ticks cannot be overlapped by the curves, which are clipped to the axis box: this emulates the
operator choosing clean ticks or tick labels, as the atlas does."""
function frame_and_ticks!(cv, xt, yt; w = 1.0, len = 2.5)
    m = hires(cv)
    stroke!(m, cv, [cv.XL, cv.XR, cv.XR, cv.XL, cv.XL], [cv.YT, cv.YT, cv.YB, cv.YB, cv.YT], w; clip = false)
    for u in xt
        X = px_of(cv, u); stroke!(m, cv, [X, X], [cv.YB, cv.YB + len], w; clip = false)
    end
    for v in yt
        Y = py_of(cv, v); stroke!(m, cv, [cv.XL, cv.XL - len], [Y, Y], w; clip = false)
    end
    composite!(cv, m, (0.0, 0.0, 0.0))
end

# ---------------------------------------------------------------- tick calibration from the image
# The atlas reads tick pixels as sub-pixel centroids of the dark tick marks. Here: darkness-weighted
# centroid along the axis of the tick pixels between the frame line and the tick tip, within +-3 px of
# the nominal position (the operator knows roughly where each tick is).
# Darkness in each row (column) of the window is taken relative to its lightest pixel, so the
# partly covered frame line does not pull the centroid.
function tick_centroid(img, along_x::Bool, nominal, frame)
    f = round(Int, frame)
    if along_x
        rows = (f + 1):(f + 3); cols = (round(Int, nominal) - 3):(round(Int, nominal) + 3)
        w = [(1 - lum(img[r, c])) for r in rows, c in cols]
        w = max.(w .- minimum(w; dims = 2), 0)
        return sum(w .* reshape(collect(cols), 1, :)) / sum(w)
    else
        cols = (f - 3):(f - 1); rows = (round(Int, nominal) - 3):(round(Int, nominal) + 3)
        w = [(1 - lum(img[r, c])) for r in rows, c in cols]
        w = max.(w .- minimum(w; dims = 1), 0)
        return sum(w .* collect(rows)) / sum(w)
    end
end

# ---------------------------------------------------------------- curve shapes (normalised display coordinates)
const U = collect(range(0.0, 1.0; length = 1501))

growth(rng; v0 = 0.05 + 0.3rand(rng), v1 = v0 + 0.3 + (0.92 - v0 - 0.3) * rand(rng), r = 0.3 + 4rand(rng)) =
    u -> v0 + (v1 - v0) * (exp(r * u) - 1) / (exp(r) - 1)
function kink(rng; v0 = 0.05 + 0.2rand(rng))     # treatment window: slow rise, fast rise to a cusp, decay
    a = 0.3 + 0.3rand(rng); b = a + 0.1 + 0.15rand(rng); h = 0.15 + (0.9 - v0 - 0.2) * rand(rng)
    τ = 0.05 + 0.2rand(rng); lvl = 0.3h
    pre(u) = v0 + 0.15 * (1 - exp(-u / 0.2))
    return u -> u <= a ? pre(u) : u <= b ? pre(a) + h * ((u - a) / (b - a))^1.5 : pre(a) + lvl + (h - lvl) * exp(-(u - b) / τ)
end
function dip_peak(rng)                             # Ouerdani P: fast dip, regrowth to a peak, slow decline
    v0 = 0.4 + 0.3rand(rng); d = 0.1 + 0.15rand(rng); τd = 0.01 + 0.04rand(rng)
    up = 0.15 + 0.2rand(rng); h = 0.1 + 0.2rand(rng); w = 0.06 + 0.1rand(rng); fin = 0.03 + 0.15rand(rng)
    return u -> fin + (v0 - fin) * exp(-u / 0.6) - d * (1 - exp(-u / τd)) * exp(-u / 0.12) + h * exp(-((u - up) / w)^2)
end
decay_from_above(rng) = (v0 = 1.3 + 0.5rand(rng); fin = 0.02 + 0.15rand(rng); τ = 0.08 + 0.2rand(rng);
                         u -> fin + (v0 - fin) * exp(-u / τ))
function spike(rng)                                # bolus: near-vertical rise, decay (steep features)
    us = 0.1 + 0.6rand(rng); vb = 0.05 + 0.3rand(rng); h = 0.3 + (0.9 - vb - 0.3) * rand(rng)
    τr = 0.001 + 0.004rand(rng); τd = 0.02 + 0.15rand(rng)
    return u -> u < us ? vb : vb + h * (1 - exp(-(u - us) / τr)) * exp(-(u - us) / τd)
end
function repeated(rng; ncyc = rand(rng, 3:8))      # Fornari progenitors: sharp drop after each dose, recovery
    vb = 0.6 + 0.3rand(rng); d = 0.3 + 0.4rand(rng); τ1 = 0.003 + 0.01rand(rng); τ2 = 0.03 + 0.05rand(rng)
    starts = range(0.03 + 0.05rand(rng), 0.85; length = ncyc)
    return u -> max(vb - sum(u > s ? d * (1 - exp(-(u - s) / τ1)) * exp(-(u - s) / τ2) : 0.0 for s in starts), 0.05)
end
function nadir(rng)                                # Fornari neutrophils: dip and recovery, gentle
    vb = 0.5 + 0.3rand(rng); d = 0.15 + 0.3rand(rng); un = 0.15 + 0.4rand(rng); w = 0.08 + 0.15rand(rng); ph = 3rand(rng)
    return u -> vb - d * exp(-((u - un) / w)^2) + 0.05 * sin(6u + ph)
end
# PK, log axis: biexponential in data units, returned in display units for an axis of `dec` decades
function pk_log(rng; top = 0.85 + 0.12rand(rng), dec = 4.0)
    α = 5 + 35rand(rng); β = 1 + 15rand(rng); f = 0.1 + 0.6rand(rng)
    return u -> top + log10(f * exp(-α * u) + (1 - f) * exp(-β * u)) / dec
end
uptake_lin(rng) = (vmax = 0.2 + 0.75rand(rng); k = 2 + 15rand(rng); u -> vmax * (1 - exp(-k * u)))

# ---------------------------------------------------------------- panel types
# geometry (axis box in px) from results/panels.csv; line widths from results/line_widths.csv
const PANELS = [
    (name = "hut_fig2_single", model = "hutchinson_2016", w = 263.9, h = 115.4, ylog = false, rule = "colour_trace"),
    (name = "hut_fig2_five", model = "hutchinson_2016", w = 263.9, h = 115.4, ylog = false, rule = "colour_trace"),
    (name = "hut_fig6", model = "hutchinson_2016", w = 238.2, h = 88.0, ylog = false, rule = "dark_magenta_small"),
    (name = "ouerdani_fig4", model = "ouerdani_2015", w = 233.0, h = 235.0, ylog = false, rule = "black_solid_dashed"),
    (name = "stroh_fig3", model = "stroh_2019", w = 296.2, h = 232.9, ylog = true, rule = "hue_track_log_markers"),
    (name = "stroh_fig4_log", model = "stroh_2019", w = 180.0, h = 149.5, ylog = true, rule = "hue_track_log_markers"),
    (name = "stroh_fig4_lin", model = "stroh_2019", w = 180.0, h = 149.0, ylog = false, rule = "hue_track_log_markers"),
    (name = "fornari_fig3", model = "fornari_2019", w = 137.3, h = 125.6, ylog = true, rule = "red_mask_log_bands"),
    (name = "fornari_fig4", model = "fornari_2019", w = 173.6, h = 155.7, ylog = true, rule = "red_mask_log_bands"),
]

# axes for a panel: x linear 0..1 data units; y linear 0..1 or log over `dec` decades from 1e-2
function axes_for(p; dec = 4.0)
    xa = Axis("x", "-", :linear, 0.0, 1.0)
    ya = p.ylog ? Axis("y", "-", :log, 1e-2, 1e-2 * 10^dec) : Axis("y", "-", :linear, 0.0, 1.0)
    return xa, ya
end
ydata(ya, v) = ya.scale === :log ? 10^(log10(ya.lo) + v * (log10(ya.hi) - log10(ya.lo))) : ya.lo + v * (ya.hi - ya.lo)

# one scene: list of curves (shape function, u range, colour, width, dash, tracer settings)
struct Curve
    f::Function
    u1::Float64; u2::Float64
    colour::NTuple{3, Float64}
    width::Float64
    dash::Union{Nothing, Tuple{Float64, Float64}}
    tag::String
    family::String
end
truth(c::Curve) = (u = filter(x -> c.u1 <= x <= c.u2, U); (u, c.f.(u)))

function scene(p, rng)
    fam = String[]
    if p.name == "hut_fig2_single"
        k = rand(rng, 1:3)
        f, nm = k == 1 ? (growth(rng), "growth") : k == 2 ? (kink(rng), "kink") : (spike(rng), "spike")
        return [Curve(f, 0.05, 0.95, MATLAB["blue"], 2.0, nothing, "blue", nm)], (;)
    elseif p.name == "hut_fig2_five"
        k = rand(rng, 1:3)
        v0 = 0.05 + 0.25rand(rng)
        cs = Curve[]
        for (i, c) in enumerate(("blue", "orange", "yellow", "purple", "green"))
            f, nm = k == 1 ? (growth(rng; v0), "growth") : k == 2 ? (kink(rng; v0), "kink") : (spike(rng), "spike")
            push!(cs, Curve(f, 0.05, 0.95, MATLAB[c], 2.0, nothing, c, nm))
        end
        return cs, (;)
    elseif p.name == "hut_fig6"
        T = rand(rng) < 0.5 ? kink(rng; v0 = 0.4 + 0.2rand(rng)) : dip_peak(rng)
        V = growth(rng; v0 = 0.3 + 0.3rand(rng), v1 = 0.4 + 0.5rand(rng), r = -2 + 4rand(rng))
        return [Curve(T, 0.15, 0.82, (0.0, 0.0, 0.0), 1.4, nothing, "black", "T"),
                Curve(V, 0.15, 0.82, (1.0, 0.0, 1.0), 1.0, nothing, "magenta", "V")], (; window = sort(0.3 .+ 0.5 .* rand(rng, 2)))
    elseif p.name == "ouerdani_fig4"
        # one crossing at most, K above P before it (the situation the Ouerdani rule is written for)
        P, K = dip_peak(rng), decay_from_above(rng)
        while true
            d = [P(u) - K(u) for u in U if K(u) < 1.0]
            count(k -> d[k] * d[k + 1] < 0, 1:length(d) - 1) <= 1 && (isempty(d) || d[1] < 0) && break
            P, K = dip_peak(rng), decay_from_above(rng)
        end
        return [Curve(P, 0.0, 1.0, (0.0, 0.0, 0.0), 1.4, nothing, "P", "dip_peak"),
                Curve(K, 0.0, 1.0, (0.0, 0.0, 0.0), 1.4, (4.0, 3.0), "K", "decay_dashed")], (;)
    elseif startswith(p.name, "stroh")
        pal = p.name == "stroh_fig3" ? (rand(rng) < 0.5 ? ["blue", "orange", "green", "lblue"] : ["blue", "yellow", "purple"]) :
              ["blue", "orange", "yellow", "purple", "green"]
        dec = p.name == "stroh_fig3" ? 4.0 : 6.0
        top = 0.85 + 0.1rand(rng)
        cs = [p.ylog ? Curve(pk_log(rng; top, dec), 0.03, 0.93, MATLAB[c], 2.0, nothing, c, "pk_log") :
                       Curve(uptake_lin(rng), 0.0, 0.95, MATLAB[c], 2.0, nothing, c, "uptake") for c in pal]
        return cs, (; markers = p.name == "stroh_fig3", dec)
    else
        f, nm = p.name == "fornari_fig3" ? (repeated(rng), "repeated") : (nadir(rng), "nadir")
        return [Curve(f, 0.0, 0.97, (0.85, 0.1, 0.1), 1.5, nothing, "red", nm)], (; bands = true)
    end
end

# ---------------------------------------------------------------- render one scene
function render(p, curves, extra, rng)
    cv = Canvas(p.w, p.h, rng)
    pix(c) = (u = filter(x -> c.u1 <= x <= c.u2, U); (px_of.(Ref(cv), u), py_of.(Ref(cv), c.f.(u))))
    if p.name == "hut_fig6"          # shaded treatment windows behind the lines
        a, b = extra.window
        m = hires(cv); fill_between!(m, cv, px_of(cv, a), px_of(cv, b), X -> cv.YT, X -> cv.YB); composite!(cv, m, (0.3, 1.0, 1.0); a = 0.7)
        m = hires(cv); fill_between!(m, cv, px_of(cv, a - 0.2), px_of(cv, a), X -> py_of(cv, 0.95), X -> cv.YB); composite!(cv, m, (1.0, 1.0, 0.4); a = 0.6)
        # drug concentration C (blue dashed, right axis) near the bottom, spiking inside the window
        u = collect(range(0.15, 0.82; length = 1500))
        cdrug = [a <= x <= b ? 0.02 + 0.9 * exp(-mod(x - a, 0.03) / 0.01) : 0.005 for x in u]
        m = hires(cv); stroke!(m, cv, px_of.(Ref(cv), u), py_of.(Ref(cv), cdrug), 1.2; dash = (3.0, 2.0)); composite!(cv, m, (0.0, 0.0, 1.0))
    end
    if get(extra, :bands, false)     # Fornari: red line first, translucent percentile bands on top, black markers
        c = curves[1]; X, Y = pix(c)
        m = hires(cv); stroke!(m, cv, X, Y, c.width); composite!(cv, m, c.colour)
        uu = filter(x -> c.u1 <= x <= c.u2, U)
        for k in 1:8
            half = 0.02 + 0.035k
            lo = c.f.(uu) .- half; hi = c.f.(uu) .+ half
            Xs = px_of.(Ref(cv), uu)
            top = X -> py_of(cv, interp(Xs, hi, X)); bot = X -> py_of(cv, interp(Xs, lo, X))
            m = hires(cv); fill_between!(m, cv, Xs[1], Xs[end], top, bot); composite!(cv, m, (0.3, 0.3, 1.0); a = 0.08)
        end
        m = hires(cv)
        for _ in 1:rand(rng, 30:80)
            u = c.u1 + (c.u2 - c.u1) * rand(rng); v = c.f(u) + 0.12 * randn(rng)
            ring!(m, cv, px_of(cv, u), py_of(cv, v), 2.2, 0.8; shape = rand(rng, (:diamond, :square, :circle)))
        end
        composite!(cv, m, (0.0, 0.0, 0.0))
    else
        phase = 10rand(rng)
        for c in curves
            X, Y = pix(c)
            m = hires(cv); stroke!(m, cv, X, Y, c.width; dash = c.dash, phase); composite!(cv, m, c.colour)
        end
    end
    if get(extra, :markers, false)   # Stroh Figure 3: hollow circles in the line colours, densest early
        for c in curves
            m = hires(cv)
            for _ in 1:rand(rng, 8:16)
                u = rand(rng) < 0.5 ? 0.04 * rand(rng) : c.u1 + (c.u2 - c.u1) * rand(rng)^1.5
                v = c.f(clamp(u, c.u1, c.u2)) + 0.03 * randn(rng)
                ring!(m, cv, px_of(cv, u), py_of(cv, v), 2.0, 0.9)
            end
            composite!(cv, m, c.colour)
        end
    end
    # four ticks per axis, evenly spaced from a random start
    tk() = (a = 0.03 + 0.12rand(rng); b = 0.85 + 0.12rand(rng); collect(range(a, b; length = 4)))
    xt = tk(); yt = tk()
    frame_and_ticks!(cv, xt, yt)
    return cv, xt, yt
end
function interp(xs, ys, x)
    x <= xs[1] && return ys[1]; x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x); t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return ys[i] + t * (ys[i + 1] - ys[i])
end

# ---------------------------------------------------------------- trace with the atlas rules
# As in the atlas: two ticks calibrate an axis and the others are checks. The operator keeps the pair
# whose calibration places the unused ticks best (smallest largest residual), so a tick spoiled by a
# curve or marker passing over it is not used.
function best_pair(t, c)
    best = (Inf, 1, 2)
    for i in 1:3, j in (i + 1):4
        pred(k) = c[i] + (t[k] - t[i]) * (c[j] - c[i]) / (t[j] - t[i])
        res = maximum(abs(pred(k) - c[k]) for k in 1:4 if k != i && k != j)
        res < best[1] && (best = (res, i, j))
    end
    return best[2], best[3]
end
function calibrations(img, cv, xt, yt, xa, ya)
    cx = [tick_centroid(img, true, px_of(cv, u), cv.YB) for u in xt]
    cy = [tick_centroid(img, false, py_of(cv, v), cv.XL) for v in yt]
    i, j = best_pair(xt, cx); k, l = best_pair(yt, cy)
    xcal = [cx[i], xt[i] * (xa.hi - xa.lo) + xa.lo, cx[j], xt[j] * (xa.hi - xa.lo) + xa.lo]
    ycal = [cy[k], ydata(ya, yt[k]), cy[l], ydata(ya, yt[l])]
    # largest error of the fitted calibration over the axis box (pixels)
    ex = maximum(abs(cx[i] + (u - xt[i]) * (cx[j] - cx[i]) / (xt[j] - xt[i]) - px_of(cv, u)) for u in (0.0, 1.0))
    ey = maximum(abs(cy[k] + (v - yt[k]) * (cy[l] - cy[k]) / (yt[l] - yt[k]) - py_of(cv, v)) for v in (0.0, 1.0))
    return xcal, ycal, max(ex, ey)
end
# exact calibration (true tick pixel positions), for the tracing-only error
exact_cal(cv, xa, ya) = ([px_of(cv, 0.0), xa.lo, px_of(cv, 1.0), xa.hi], [py_of(cv, 0.0), ydata(ya, 0.0), py_of(cv, 1.0), ydata(ya, 1.0)])

# box inside the frame: `m` pixels from the frame lines
inner_box(cv, m) = [ceil(Int, cv.YT + m), floor(Int, cv.YB - m), ceil(Int, cv.XL + m), floor(Int, cv.XR - m)]

# Ouerdani rule (trace_fig4.jl): dark runs (lum < 0.6) per column, darkness-weighted centroid of each run;
# before the crossing column P = lowest run, K = the others; after it P = topmost run, K = the others;
# the crossing column and its neighbours are skipped.
function trace_ouerdani(img, box, cross_col)
    r1, r2, c1, c2 = box
    P = Tuple{Int, Float64}[]; K = Tuple{Int, Float64}[]
    for c in c1:c2
        rs = [r for r in r1:r2 if lum(img[r, c]) < 0.6]
        isempty(rs) && continue
        cen = [sum(r * (1 - lum(img[r, c])) for r in a:b) / sum(1 - lum(img[r, c]) for r in a:b) for (a, b) in runs(rs)]
        abs(c - cross_col) <= 1 && continue
        if c < cross_col
            push!(P, (c, cen[end])); for x in cen[1:end-1]; push!(K, (c, x)); end
        else
            push!(P, (c, cen[1])); for x in cen[2:end]; push!(K, (c, x)); end
        end
    end
    return P, K
end
thin(v, n) = v[unique(round.(Int, range(1, length(v); length = min(n, length(v)))))]

# Fornari rule (digitise_run.jl trace_red): red > 0.4 and red - green > 0.20, median row per column
function trace_fornari(img, box)
    r1, r2, c1, c2 = box
    px = Float64[]; py = Float64[]
    for c in c1:c2
        rs = [r for r in r1:r2 if red(img[r, c]) > 0.4 && red(img[r, c]) - green(img[r, c]) > 0.20]
        isempty(rs) && continue
        push!(px, c); push!(py, median(rs))
    end
    return px, py
end

function digitise(p, curves, extra, cv, path, xcal, ycal, xa, ya)
    img = loadimg(path)
    ysc = ya.scale === :log ? "log" : "linear"
    fx = calib(xcal..., "linear"); fy = calib(ycal..., ysc)
    out = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    root, file = dirname(path), basename(path)
    if startswith(p.name, "hut_fig2")
        box = inner_box(cv, 1.5)
        for c in curves
            spec = Dict{String, Any}("id" => c.tag, "image" => file, "box" => box, "colour" => HUT_TARGETS[c.tag], "tol" => 0.2,
                                     "xcal" => xcal, "ycal" => ycal, "n" => 60)
            try
                x, y = trace_rects(spec, root)[1:2]; out[c.tag] = (x, y)
            catch e
                e isa ErrorException || rethrow()
            end
        end
    elseif p.name == "hut_fig6"
        box = inner_box(cv, 3)
        pr = [ceil(Int, px_of(cv, curves[1].u1)), floor(Int, px_of(cv, curves[1].u2))]
        specs = [("black", Dict{String, Any}("colour" => "dark", "tol" => 0.2, "dark_thr" => 0.6, "outlier_px" => 4.0)),
                 ("magenta", Dict{String, Any}("colour" => [0.93, 0.45, 0.93], "tol" => 0.25, "outlier_px" => 3.0))]
        for (tag, s) in specs
            spec = merge(s, Dict{String, Any}("id" => tag, "image" => file, "box" => box, "xcal" => xcal, "ycal" => ycal, "n" => 70,
                                              "px_range" => pr, "outlier_window" => 20))
            try
                x, y = trace_rects(spec, root)[1:2]; out[tag] = (x, y)
            catch e
                e isa ErrorException || rethrow()
            end
        end
    elseif p.name == "ouerdani_fig4"
        P, K = curves
        # crossing column (P = K): known to the operator from the image, here from the true curves
        uu = U[2:end]; d = P.f.(uu) .- K.f.(uu); i = findfirst(k -> d[k] * d[k + 1] <= 0 && K.f(uu[k]) < 1, 1:length(d) - 1)
        cross = i === nothing ? 10^6 : round(Int, px_of(cv, uu[i]))
        Pp, Kp = trace_ouerdani(img, inner_box(cv, 3), cross)
        for (tag, pts, n) in (("P", Pp, 90), ("K", Kp, 70))
            isempty(pts) && continue
            pts = thin(sort(pts), n)
            out[tag] = (fx.(first.(pts)), fy.(last.(pts)))
        end
    elseif startswith(p.name, "stroh")
        box = inner_box(cv, 1.5)
        fourcol = length(curves) == 4
        for c in curves
            spec = Dict{String, Any}("id" => c.tag, "image" => file, "box" => box, "colour" => collect(c.colour), "tol" => 0.4,
                                     "xcal" => xcal, "ycal" => ycal, "xscale" => "linear", "yscale" => ysc, "n" => 60)
            if fourcol && c.tag in ("blue", "lblue")
                # tight, non-overlapping hue windows as in the atlas (209 and 197 there); here centred so that the
                # drawn hues (204 and 198) fall inside their own window
                spec["hue"] = c.tag == "blue" ? 207.0 : 195.0; spec["hue_tol"] = 5.0
            end
            try
                px, py = robust_trace(spec, root)
                keep = get(extra, :markers, false) ? fx.(px) .>= 1 / 31 : trues(length(px))   # min_x = 1 d on a 31 d axis (Figure 3)
                px, py = px[keep], py[keep]
                isempty(px) && continue
                idx = unique(round.(Int, range(1, length(px); length = min(60, length(px)))))
                out[c.tag] = (fx.(px[idx]), fy.(py[idx]))
            catch e
                e isa ErrorException || rethrow()
            end
        end
    else
        px, py = trace_fornari(img, inner_box(cv, 1.5))
        if !isempty(px)
            idx = unique(round.(Int, range(1, length(px); length = min(140, length(px)))))
            out["red"] = (fx.(px[idx]), fy.(py[idx]))
        end
    end
    return out
end

# ---------------------------------------------------------------- driver
function run_one(p, rep, seed, tmpdir; keep = false)
    rng = MersenneTwister(hash((seed, p.name, rep)))
    curves, extra = scene(p, rng)
    dec = haskey(extra, :dec) ? extra.dec : (p.name == "fornari_fig3" ? 3.0 : p.name == "fornari_fig4" ? 3.0 : 4.0)
    xa, ya = axes_for(p; dec)
    cv, xt, yt = render(p, curves, extra, rng)
    path = joinpath(tmpdir, "$(p.name)_$(rep).jpg")
    grey = p.name == "ouerdani_fig4"
    save(path, grey ? Gray.(cv.img) : cv.img; quality = 92)
    xcal, ycal, calerr = calibrations(loadimg(path), cv, xt, yt, xa, ya)
    traced = digitise(p, curves, extra, cv, path, xcal, ycal, xa, ya)
    traced0 = digitise(p, curves, extra, cv, path, exact_cal(cv, xa, ya)..., xa, ya)
    rows = NamedTuple[]
    nan = (rms = NaN, p95 = NaN, max = NaN)
    sc(t, c, tx, ty) = haskey(t, c.tag) && length(t[c.tag][1]) >= 3 ? score(t[c.tag]..., tx, ty, xa, ya; id = c.tag) : nan
    for c in curves
        u, v = truth(c)
        tx = u .* (xa.hi - xa.lo) .+ xa.lo; ty = ydata.(Ref(ya), v)
        r = sc(traced, c, tx, ty); r0 = sc(traced0, c, tx, ty)
        n = haskey(traced, c.tag) ? length(traced[c.tag][1]) : 0
        push!(rows, (; panel = p.name, model = p.model, rule = p.rule, rep, curve = c.tag, family = c.family, n,
                     rms = r.rms, p95 = r.p95, max = r.max, rms_trace = r0.rms, p95_trace = r0.p95, traced = n >= 3, calerr_px = calerr))
    end
    keep || rm(path)
    return rows
end

function main(reps, seed)
    outdir = joinpath(@__DIR__, "results"); mkpath(outdir)
    exdir = joinpath(@__DIR__, "examples"); mkpath(exdir)
    tmp = mktempdir()
    jobs = [(p, r) for p in PANELS for r in 1:reps]
    results = Vector{Vector{NamedTuple}}(undef, length(jobs))
    Threads.@threads for k in eachindex(jobs)
        p, r = jobs[k]
        results[k] = run_one(p, r, seed, tmp; keep = r <= 2)
    end
    for p in PANELS, r in 1:2       # keep two synthetic examples per panel type
        f = joinpath(tmp, "$(p.name)_$(r).jpg"); isfile(f) && cp(f, joinpath(exdir, "$(p.name)_$(r).jpg"); force = true)
    end
    rows = reduce(vcat, results)
    open(joinpath(outdir, "synthetic_runs.csv"), "w") do io
        println(io, "panel,model,rule,rep,curve,family,n,rms,p95,max,rms_trace_only,p95_trace_only,traced,calerr_px")
        for r in rows
            @printf(io, "%s,%s,%s,%d,%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.3f\n", r.panel, r.model, r.rule, r.rep, r.curve, r.family, r.n, r.rms, r.p95, r.max, r.rms_trace, r.p95_trace, r.traced, r.calerr_px)
        end
    end
    println(length(rows), " traced curves written")
end

if abspath(PROGRAM_FILE) == @__FILE__
    reps = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200
    seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2026
    main(reps, seed)
end
