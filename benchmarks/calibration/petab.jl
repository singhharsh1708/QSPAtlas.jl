# PEtab.jl maximum-likelihood recalibration of hutchinson_2016 against Figure 1d (group means).
# Run from QSPAtlas.jl/:  julia --startup-file=no --project=benchmarks/calibration benchmarks/calibration/petab.jl
using Dates
println("start ", now())
include(joinpath(@__DIR__, "common.jl"))
using PEtab, Fides, Optim, Distributions, Random, StatsBase, LinearAlgebra, TOML, Printf
import Distributions: LogNormal
using ForwardDiff

NSTART = parse(Int, get(ENV, "NSTART", "50"))
T0MODE = get(ENV, "T0MODE", "fixed")     # fixed | shared | group  (sensitivity runs)

# NOTE: PEtab.jl works in time since day 38 (tau = t - 38): SBMLImporter's time-event initialisation adds a tstop at
# min(event times)/2, which is only valid when the simulation starts at 0 (otherwise 'tstop behind the current time').
# The ODE has no explicit time dependence outside the window events, so the shift changes nothing.
# ---- observable, conditions, events, parameters --------------------------------------------------
@parameters sigma delta_V_bev delta_V_van Nmax T0_ctrl T0_bev T0_van T0_all
obs = PEtabObservable(:obs_T, :obsT, sigma; distribution = LogNormal)

t0sym(g) = T0MODE == "fixed" ? T0_FIXED : T0MODE == "shared" ? T0_all : (g == "Control" ? T0_ctrl : g == "Bevacizumab" ? T0_bev : T0_van)
conds = [
    PEtabCondition(:ctrl, :delta_V => 0.0, :Nwin => 1.0, :T0 => t0sym("Control"); t0 = 0.0),
    PEtabCondition(:bev, :delta_V => delta_V_bev, :Nwin => Nmax, :T0 => t0sym("Bevacizumab"); t0 = 0.0),
    PEtabCondition(:van, :delta_V => delta_V_van, :Nwin => Nmax, :T0 => t0sym("Vanucizumab"); t0 = 0.0),
]
tt = t
events = [
    PEtabEvent(tt == PUB.tnorm1 - T_START, :Nw => Nwin),   # window opens: N = Nmax (condition-specific Nwin, 1 for control)
    PEtabEvent(tt == PUB.tnorm2 - T_START, :Nw => 1.0),    # window closes
]

pest = [
    PEtabParameter(:alpha_T; lb = 1.0e-3, ub = 5.0, value = PUB.alpha_T),
    PEtabParameter(:alpha_V; lb = 1.0e-3, ub = 5.0, value = PUB.alpha_V),
    PEtabParameter(:K; lb = 0.05, ub = 50.0, value = PUB.K),
    PEtabParameter(:Nmax; lb = 1.0, ub = 1.0e4, value = PUB.Nmax),
    PEtabParameter(:delta_V_bev; lb = 1.0e-3, ub = 3.0, value = PUB.delta_V_bev),
    PEtabParameter(:delta_V_van; lb = 1.0e-3, ub = 3.0, value = PUB.delta_V_van),
    PEtabParameter(:sigma; lb = 1.0e-3, ub = 3.0, value = 0.2),
]
if T0MODE == "shared"
    push!(pest, PEtabParameter(:T0_all; lb = 10.0, ub = 300.0, value = 70.0))
elseif T0MODE == "group"
    append!(pest, [PEtabParameter(:T0_ctrl; lb = 10.0, ub = 300.0, value = 70.0), PEtabParameter(:T0_bev; lb = 10.0, ub = 300.0, value = 70.0), PEtabParameter(:T0_van; lb = 10.0, ub = 300.0, value = 70.0)])
end

df = load_data()
simid = Dict("Control" => "ctrl", "Bevacizumab" => "bev", "Vanucizumab" => "van")
meas = DataFrame(simulation_id = [simid[g] for g in df.group], obs_id = "obs_T", time = Float64.(df.time) .- T_START, measurement = df.volume)

model = PEtabModel(FITSYS, [obs], meas, pest; simulation_conditions = conds, events = events)
prob = PEtabODEProblem(model; odesolver = ODESolver(Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10))
println(describe(prob))

x0 = get_x(prob)
println("nominal x ", x0)
println("nllh at nominal (published values, sigma = 0.2) ", prob.nllh(x0))

# ---- published parameter vector with its best sigma ------------------------------------------------
function sigma_best(prob, xpub)
    f(ls) = (x = copy(xpub); x.log10_sigma = ls; prob.nllh(x))
    r = Optim.optimize(f, log10(1.0e-3), log10(3.0), Optim.Brent())
    x = copy(xpub); x.log10_sigma = Optim.minimizer(r)
    return x, Optim.minimum(r)
end
xpub = copy(x0)
xpub, nllh_pub = sigma_best(prob, xpub)
println("published parameters, sigma re-optimised: nllh = ", nllh_pub, " sigma = ", 10^xpub.log10_sigma)

# ---- multistart ---------------------------------------------------------------------------------------
rng = MersenneTwister(20261104)
if get(ENV, "PROFILES_ONLY", "0") == "1"
    m = TOML.parsefile(joinpath(@__DIR__, "results", "petab_mle.toml"))
    xbest = copy(x0)
    for nm in keys(xbest)
        xbest[nm] = log10(m[String(nm)[7:end]])
    end
    fbest = prob.nllh(xbest)
    ms = nothing
else
t_ms = @elapsed ms = calibrate_multistart(rng, prob, Fides.BFGS(), NSTART)
println("multistart done in ", round(t_ms; digits = 1), " s")
println(ms)
xbest = ms.xmin
fbest = ms.fmin
end
println("best nllh ", fbest, "  x ", xbest)

if ms !== nothing
    # collect runs
    runs = DataFrame(run = 1:length(ms.runs), fmin = [r.fmin for r in ms.runs], niter = [r.niterations for r in ms.runs])
    for nm in keys(xbest)
        runs[!, nm] = [r.xmin[nm] for r in ms.runs]
    end
    sort!(runs, :fmin)
    CSV.write(joinpath(@__DIR__, "results", "petab_multistart$(T0MODE == "fixed" ? "" : "_T0" * T0MODE).csv"), runs)
    println(first(runs, 12))
    println("runs within 0.01 of best: ", count(<(fbest + 0.01), runs.fmin), " / ", nrow(runs),
        "; within 2 (LR 4 units): ", count(<(fbest + 2), runs.fmin))
end

# ---- Hessian-based covariance (log10 scale) -------------------------------------------------------------
H = prob.hess(xbest)
ev = eigvals(Symmetric(H))
println("Hessian eigenvalues (log10 scale) ", ev)
cov = try inv(Symmetric(H)) catch; nothing end
if cov !== nothing
    sd = sqrt.(abs.(diag(cov)))
    println("param names ", keys(xbest))
    corr = cov ./ (sd * sd')
    println("sd log10 ", Dict(zip(keys(xbest), sd)))
    println("corr (log10 scale):"); show(stdout, "text/plain", round.(corr; digits = 2)); println()
end

# ---- profile likelihood (own implementation: fix one parameter on a grid, re-optimise the rest) ----------
function profile_param(prob, xbest, nm; npts = 25, span = 1.0, lbs = nothing)
    # grid in log10 around the MLE; at each grid value the parameter is fixed with remake(prob; parameters) and the
    # other parameters are re-optimised with Fides (exact ForwardDiff Hessian), warm-started from the neighbouring point
    i = findfirst(==(nm), collect(keys(xbest)))
    lb = collect(prob.lower_bounds)[i]; ub = collect(prob.upper_bounds)[i]
    xb = xbest[nm]
    grid_lo = max(lb, xb - span); grid_hi = min(ub, xb + span)
    grid = sort(unique(vcat(range(grid_lo, xb; length = npts ÷ 2 + 1), range(xb, grid_hi; length = npts ÷ 2 + 1))))
    pname = Symbol(String(nm)[7:end])
    res = Tuple{Float64, Float64}[]
    for dir in (1, -1)
        gs = dir == 1 ? filter(>(xb), grid) : reverse(filter(<(xb), grid))
        xw = copy(xbest)
        for gv in gs
            pr = remake(prob; parameters = [pname => 10.0^gv])
            x0 = get_x(pr)
            for k in keys(x0)
                x0[k] = xw[k]
            end
            r = calibrate(pr, x0, Fides.CustomHessian())
            for k in keys(x0)
                xw[k] = r.xmin[k]
            end
            push!(res, (gv, r.fmin))
        end
    end
    push!(res, (xb, prob.nllh(xbest)))
    return sort(unique(res); by = first)
end
using ComponentArrays
prof = Dict{Symbol, Vector{Tuple{Float64, Float64}}}()
if get(ENV, "PROFILES", "1") == "1" && T0MODE == "fixed"
    for nm in keys(xbest)
        t_p = @elapsed prof[nm] = profile_param(prob, xbest, nm; span = nm == :log10_Nmax ? 4.0 : 1.0, npts = nm == :log10_Nmax ? 31 : 25)
        println("profile ", nm, " done in ", round(t_p; digits = 1), " s")
    end
    prows = DataFrame(parameter = String[], log10_value = Float64[], value = Float64[], nllh = Float64[])
    for (nm, v) in prof, (gv, fv) in v
        push!(prows, (string(nm), gv, 10^gv, fv))
    end
    CSV.write(joinpath(@__DIR__, "results", "petab_profiles.csv"), prows)
end

# ---- predictions on a grid and per-observation residuals ------------------------------------------------
prob_sim = make_prob()
function theta(x)
    (alpha_T = 10^x.log10_alpha_T, alpha_V = 10^x.log10_alpha_V, K = 10^x.log10_K, Nmax = 10^x.log10_Nmax,
        delta_V_bev = 10^x.log10_delta_V_bev, delta_V_van = 10^x.log10_delta_V_van, sigma = 10^x.log10_sigma)
end
θbest = theta(xbest); θpub_fit = theta(xpub)
# cross-check PEtab simulated values against the plain remake simulation used by Turing
sv = prob.simulated_values(xbest)
chk = Float64[]
for g in GROUPS
    idx = findall(==(g), df.group)
    ys = simulate_group(prob_sim, g, θbest, Float64.(df.time[idx]))
    append!(chk, ys)
end
# PEtab orders measurements as in `meas`; df rows are ordered by group blocks, so compare directly
println("max |PEtab sim - remake sim| = ", maximum(abs.(sv .- chk)))

pred = DataFrame(group = String[], time = Float64[], published = Float64[], mle = Float64[])
for g in GROUPS
    yp = simulate_group(prob_sim, g, (; (k => PUB[k] for k in PARNAMES[1:6])...), TGRID)
    ym = simulate_group(prob_sim, g, θbest, TGRID)
    for (a, b, c) in zip(TGRID, yp, ym)
        push!(pred, (g, a, b, c))
    end
end
CSV.write(joinpath(@__DIR__, "results", "petab_predictions$(T0MODE == "fixed" ? "" : "_T0" * T0MODE).csv"), pred)

if ms !== nothing
    est = Dict("nllh" => fbest, "nllh_published_sigma_opt" => nllh_pub, "n_obs" => nrow(df), "n_starts" => NSTART,
        "T0mode" => T0MODE,
        "runs_within_0.01" => ms === nothing ? -1 : count(<(fbest + 0.01), runs.fmin))
    for nm in keys(xbest)
        est[String(nm)[7:end]] = 10^xbest[nm]
    end
    open(joinpath(@__DIR__, "results", "petab_mle$(T0MODE == "fixed" ? "" : "_T0" * T0MODE).toml"), "w") do io
        TOML.print(io, est)
    end
end
println("LR statistic published vs MLE: ", 2 * (nllh_pub - fbest))
println("done ", now())
