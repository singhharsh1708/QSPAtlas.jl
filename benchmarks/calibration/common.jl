# Shared definitions for the hutchinson_2016 recalibration (PEtab.jl and Turing.jl scripts).
# Model: Eqs 1-5 of Hutchinson 2016 with the K typo fix (V0 = K*T0). The window factor N(t) is a parameter `Nw`
# switched by two time events (PEtab events in PEtab.jl, three ODE segments in Turing.jl), so no discontinuous
# right-hand side has to be integrated through.
using ModelingToolkit, OrdinaryDiffEq, DataFrames, CSV
using ModelingToolkit: t_nounits as t, D_nounits as D

const HERE = @__DIR__
const T_START = 38.0
const T_END = 71.0
const GROUPS = ["Control", "Bevacizumab", "Vanucizumab"]

# Table 1 final model (population values), PMC5192995
const PUB = (alpha_T = 0.109, alpha_V = 0.119, K = 1.14, Nmax = 6.7, delta_V_bev = 0.113, delta_V_van = 0.115,
    tnorm1 = 53.4, tnorm2 = 59.2, sigma = 0.185)
const T0_FIXED = 70.0      # Methods: mean size 70 mm^3 at day 38; first digitised means 70, 70, 70.9
const PARNAMES = [:alpha_T, :alpha_V, :K, :Nmax, :delta_V_bev, :delta_V_van, :sigma]

function load_data(; include_t0 = false)
    df = CSV.read(joinpath(HERE, "data", "hutchinson_fig1d.csv"), DataFrame; comment = "#")
    return include_t0 ? df : filter(r -> r.time > T_START, df)
end

@parameters alpha_T = PUB.alpha_T alpha_V = PUB.alpha_V K = PUB.K T0 = T0_FIXED delta_V = 0.0 Nwin = 1.0 Nw = 1.0 tn1 = PUB.tnorm1 tn2 = PUB.tnorm2
@variables T(t) = T0 V(t) = K * T0 obsT(t) obsNwin(t)
eqs = [
    D(T) ~ alpha_T * T * (1 - T / (Nw * V)),
    D(V) ~ alpha_V * T^(2 / 3) * V^(1 / 3) - delta_V * V,
    obsT ~ T,
    obsNwin ~ Nwin,    # keeps the window-height parameter (set per group by PEtab conditions) in the compiled system
]
@named hsys = System(eqs, t)
const FITSYS = mtkcompile(hsys)

# group settings: vessel kill parameter and whether the window acts
group_params(g, θ) = g == "Control" ? (delta_V = 0.0, Nmax = 1.0) :
    g == "Bevacizumab" ? (delta_V = θ.delta_V_bev, Nmax = θ.Nmax) : (delta_V = θ.delta_V_van, Nmax = θ.Nmax)

"""
    simulate_group(prob, g, θ, ts; T0val = T0_FIXED, tn1, tn2, solver, kw...) -> Vector (T at ts)

Simulate T(ts) for group `g` with parameters `θ` (NamedTuple with alpha_T, alpha_V, K, Nmax, delta_V_bev, delta_V_van)
as three segments (before / inside / after the window) of the MTK problem `prob` via `remake`. Works with
ForwardDiff dual numbers in `θ`.
"""
function simulate_group(prob, g, θ, ts; T0val = T0_FIXED, tn1 = PUB.tnorm1, tn2 = PUB.tnorm2,
        solver = Tsit5(), abstol = 1.0e-9, reltol = 1.0e-9)
    gp = group_params(g, θ)
    RT = promote_type(typeof(θ.alpha_T), typeof(θ.K), typeof(gp.delta_V), typeof(gp.Nmax))
    u = RT[T0val, θ.K * T0val]
    out = zeros(RT, length(ts))
    bounds = (T_START, tn1, tn2, T_END)
    Nws = (one(RT), gp.Nmax, one(RT))
    for s in 1:3
        a, b = bounds[s], bounds[s + 1]
        p = [alpha_T => θ.alpha_T, alpha_V => θ.alpha_V, delta_V => gp.delta_V, Nw => Nws[s]]
        pr = remake(prob; u0 = [T => u[1], V => u[2]], p = p, tspan = (a, b))
        sol = OrdinaryDiffEq.solve(pr, solver; abstol, reltol, saveat = Float64[], dense = true)
        for (i, tt) in enumerate(ts)
            if (tt > a || (s == 1 && tt >= a)) && tt <= b
                out[i] = sol(tt; idxs = T)
            end
        end
        u = RT[sol(b; idxs = T), sol(b; idxs = V)]
    end
    return out
end

"ODE problem for simulate_group (no initialisation problem: the segment start values are passed explicitly)."
make_prob() = ODEProblem(FITSYS, [], (T_START, T_END); build_initializeprob = false)

# Dense simulation grid for plots
const TGRID = collect(range(T_START, T_END; length = 200))
