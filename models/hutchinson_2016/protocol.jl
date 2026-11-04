# Protocol for the Hutchinson 2016 reproduction: scenarios, targets, tolerances, fixes, benchmark setups.
#
# Targets were digitised and frozen on 2026-11-04 at about 09:48 IST (reference/*.csv, reference/digitise.toml),
# before any simulation was compared with them. Tolerances below are the atlas defaults, fixed at the same time.
using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

const TOLERANCE = Tolerance()                         # RMS 0.015, P95 0.03 of the axis box (default)
const TOLERANCES = Dict{String, Tolerance}()          # no per-curve overrides

const T_START = 38.0     # day of first measurement and start of treatment (Methods, Preclinical data)
const T_END = 71.0       # last time drawn in Fig 2 and Fig 6 (Fig 4 / Fig 6 axes end at day 71)

# Fixes (only with evidence; see report.md). Ids are used by simulate_curve.
const FIXES = [
    Fix("K_is_V_over_T", :typo,
        "Initial carrying capacity V0 = K*T0 instead of the printed K = T0/V0 (i.e. V0 = T0/K).";
        evidence = "Methods before Eq. 6 prints K = T0/V0; Table 1 describes K as 'Initial ratio of V/T'; supplement MLXTRAN code (s004) has V_0 = K*X, T_0 = X; " *
                   "Fig 2b starts at V(38) = 140 = K*T0 with K = 2, T0 = 70 as stated in the caption."),
    # The next two apply to Figure 6 only (scenarios with fig = :fig6); Figure 2 is simulated exactly as before.
    Fix("resistance_clock_from_tCon", :ambiguity,
        "Figure 6: the resistance factor of Eq. 7 is exp(-lambda (t - tCon)), timed from the start of chemotherapy, instead of exp(-lambda t) with t read as days since inoculation.";
        evidence = "Eq. 7 writes exp(-lambda t) and nowhere defines the origin of t; the text says resistance builds up during chemotherapy and is likely to matter more in Fig 6(d) (long course) than in Fig 6(c). " *
                   "Both readings are textually possible; the one used is stated here. No parameter value is changed.",
        figures = ["fig6"]),
    Fix("T0_fig6_from_figure", :from_target,
        "Figure 6: T0 = 80 mm^3 (Methods mean of 70 mm^3 not used) because the curves of all four panels start at T of about 80.";
        evidence = "The paper gives no T0 for Figure 6. Value read from the starting point of the T curves in Figure 6 itself (a protocol input, not a model parameter); " *
                   "an independent rebuild from the paper alone read the same 80. Circular for this one value, hence the dagger.",
        figures = ["fig6"]),
]

# ---- scenarios ---------------------------------------------------------------------------------
fig2_base = (; alpha_T = 0.1, alpha_V = 0.09, K = 2.0, T0 = 70.0)                  # Fig 2 caption
fig6_base = (; alpha_T = 0.109, alpha_V = 0.119, K = 1.14, T0 = 70.0,              # Table 1 final model (population); T0 = 70: Methods mean size at day 38
    tnorm1 = 53.4, tnorm2 = 59.2, Nmax = 6.7, delta_V = 0.115,                      # Table 1, vanucizumab (delta_V van = 0.115)
    delta_T = 0.12, lambda = 0.08, alphaC_on = 1.0, k_C = 0.9,                      # chemotherapy section text
    fig = :fig6)                                                                    # marks scenarios the Figure 6 fixes apply to

const DV = Dict("d04" => 0.04, "d06" => 0.06, "d08" => 0.08, "d10" => 0.10, "d12" => 0.12)  # Fig 2c,d legend
const NM = Dict("N2" => 2.0, "N3" => 3.0, "N6" => 6.0, "N11" => 11.0, "N21" => 21.0)        # Fig 2e,f legend
const CHEMO = Dict("a" => (0.0, 0.0), "b" => (42.0, 49.0), "c" => (54.0, 61.0), "d" => (38.0, 59.0))  # (tCon, tCoff)

"""
    scenario(id) -> (params::NamedTuple, output::Symbol)

Parameters and plotted state (`:T` or `:V`) for reference curve `id`.
"""
function scenario(id::AbstractString)
    parts = split(id, '_')
    fig, panel = parts[1], parts[2]
    if fig == "fig2"
        out = panel in ("a", "c", "e") ? :T : :V
        p = if panel in ("a", "b")
            merge(fig2_base, (; delta_V = 0.0, Nmax = 1.0))
        elseif panel in ("c", "d")
            merge(fig2_base, (; delta_V = DV[parts[4]], Nmax = 1.0))
        else
            merge(fig2_base, (; delta_V = 0.1, Nmax = NM[parts[4]], tnorm1 = 52.0, tnorm2 = 62.0))
        end
        return p, out
    elseif fig == "fig6"
        on, off = CHEMO[panel]
        # panel a has no chemotherapy: alphaC_on stays 1 but tCon = tCoff = 0 so the drug is never on
        return merge(fig6_base, (; tCon = on, tCoff = off)), Symbol(parts[3])
    end
    error("unknown reference id $id")
end

const SYS = hutchinson_system()

"""
    breakpoint_times(p)

Times in (T_START, T_END) where the right-hand side is discontinuous: window edges, chemotherapy on/off, and
every integer day while the drug is on (the sawtooth of Eq. 9). The integrator is made to step exactly onto them.
"""
function breakpoint_times(p)
    ts = Float64[]
    get(p, :Nmax, 1.0) != 1.0 && append!(ts, [p.tnorm1, p.tnorm2])
    if get(p, :tCoff, 0.0) > get(p, :tCon, 0.0)
        append!(ts, [p.tCon, p.tCoff])
        append!(ts, collect(ceil(p.tCon):floor(p.tCoff)))
    end
    return sort(unique(filter(x -> T_START < x < T_END, ts)))
end

"""
    hutchinson_problem(p; fixes) -> (prob, callback)

ODE problem for scenario parameters `p`. The callback adds 0 to `T` at the breakpoints, which only forces
the solver to stop (and restart its derivative) on the discontinuities.
"""
function hutchinson_problem(p; fixes = Fix[])
    ids = Set(f.id for f in fixes)
    is6 = get(p, :fig, :fig2) === :fig6
    "T0_fig6_from_figure" in ids && is6 && (p = merge(p, (; T0 = 80.0)))
    "resistance_clock_from_tCon" in ids && is6 && (p = merge(p, (; t_res0 = p.tCon)))
    T0, V0 = initial_state(p.T0, p.K; K_is_V_over_T = "K_is_V_over_T" in ids)
    kv = Pair[SYS.T => T0, SYS.V => V0]
    names = Dict(:alpha_T => SYS.alpha_T, :alpha_V => SYS.alpha_V, :delta_V => SYS.delta_V, :Nmax => SYS.Nmax,
        :tnorm1 => SYS.tnorm1, :tnorm2 => SYS.tnorm2, :delta_T => SYS.delta_T, :lambda => SYS.lambda,
        :alphaC_on => SYS.alphaC_on, :k_C => SYS.k_C, :tCon => SYS.tCon, :tCoff => SYS.tCoff,
        :t_res0 => SYS.t_res0, :Noff => SYS.Noff)
    for (k, par) in names
        haskey(p, k) && push!(kv, par => Float64(p[k]))
    end
    bt = breakpoint_times(p)
    cb = isempty(bt) ? nothing : bolus_callback(SYS, SYS.T, bt, 0.0)
    prob = ODEProblem(SYS, kv, (T_START, T_END))
    return prob, cb
end

const GRID = collect(T_START:0.05:T_END)

"""
    simulate_scenario(p; fixes, grid) -> NamedTuple (t, T, V, N, C)

Solve scenario parameters `p` and sample all states and observed variables on `grid`.
"""
function simulate_scenario(p; fixes = Fix[], grid = GRID, overrides = (;))
    p = merge(p, overrides)
    prob, cb = hutchinson_problem(p; fixes)
    kw = cb === nothing ? (;) : (; callback = cb)
    sol = solve(prob, Tsit5(); abstol = 1e-9, reltol = 1e-9, kw...)
    get(v) = collect(Float64, sol(grid; idxs = v).u)
    return (; t = grid, T = get(SYS.T), V = get(SYS.V), N = get(SYS.N), C = get(SYS.C))
end

"""
    simulate_curve(ref; fixes) -> (x, y)

Simulated curve for reference `ref` on a 0.05-day grid.
"""
function simulate_curve(ref; fixes = Fix[], overrides = (;))
    p, out = scenario(ref.id)
    s = simulate_scenario(p; fixes, overrides)
    return s.t, out === :T ? s.T : s.V
end

"""
    benchmark_setups() -> Vector{BenchmarkSetup}

Problems used for the figures (with the documented fixes applied): control (Fig 2a), window with the largest
enhancement (Fig 2e, Nmax = 21) and the two chemotherapy regimens with a sawtooth drug input (Fig 6c, 6d).
"""
function benchmark_setups()
    out = BenchmarkSetup[]
    for (name, id) in [("fig2a control", "fig2_a_T"), ("fig2e Nmax=21 window", "fig2_e_T_N21"),
                       ("fig6c chemo in window", "fig6_c_T"), ("fig6d chemo 38-59", "fig6_d_T")]
        p, _ = scenario(id)
        prob, cb = hutchinson_problem(p; fixes = FIXES)
        push!(out, BenchmarkSetup(prob; name, callback = cb, saveat = collect(T_START:0.5:T_END)))
    end
    return out
end
