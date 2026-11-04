# Population pilot on an existing atlas model: Ouerdani et al. 2015, preclinical (CAKI-2
# mice) tumour growth / angiogenesis inhibition model, 100 mg/kg. 500 virtual mice.
#
# Run from QSPAtlas.jl:  julia --project=. examples/population_ouerdani.jl
#
# Inter-individual variability (IIV), log-normal, uncorrelated, taken from
#   Table 1 of the paper (preclinical IIV column, "square root of the variance (omega)", %):
#     lambda 53 %, K0 36 %, gamma 19 %, alpha 24 %, delta 42 %; b fixed to 0 (no IIV).
#   Cross-check: NONMEM $OMEGA in the supplement (PSP4-4-660-s004.txt) lists the variances
#     0.273 (lambda), 0.127 (K0), 0.0418 (gamma), 0.0506 (alpha), 0.167 (delta); their square
#     roots (0.52, 0.36, 0.20, 0.22, 0.41) agree with Table 1 within its rounding and the
#     difference between initial and final estimates. Table 1 is used.
#   The NONMEM code applies the gamma, alpha and delta etas to the whole rate (the
#   AUC-dependent gamma = gamma0 AUC^beta_gamma exp(eta)), which is the same as scaling
#   gamma0, alpha0 and delta here.
# Residual error: combined, 14 % proportional and 3 mm^3 additive (Table 1, as sd).
# Not modelled: IIV on the initial volume. The paper sets P0 to each animal's observed
#   baseline; the atlas protocol uses the typical P0 = 150 mm^3, so every mouse starts there.
using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using Statistics: quantile
import CairoMakie

Base.include(@__MODULE__, joinpath(models_dir(), "ouerdani_2015", "model.jl"))     # TUMOUR, PRECLINICAL, MOUSE_AUC
Base.include(@__MODULE__, joinpath(models_dir(), "ouerdani_2015", "protocol.jl"))  # FIXES, problem_for, ...

const TOL = (abstol = 1e-10, reltol = 1e-10)
ω = (λ = 0.53, K0 = 0.36, γ = 0.19, α = 0.24, δ = 0.42)          # Table 1 IIV as sd of eta
ηnames = [:λ, :K0, :γ, :α, :δ]
Ω = [ω[n]^2 * (i == j) for (i, n) in enumerate(ηnames), (j, _) in enumerate(ηnames)]

par = PRECLINICAL
model = PopulationModel(TUMOUR;
    routes = Route[],                       # no dosing events: treatment acts through the exposure AUC
    θ = (par..., P0 = FIXED_P0.mouse),
    ηnames, Ω,
    params = (θ, η, c) -> [λ => θ.λ * exp(η.λ), b => θ.b, γ0 => θ.γ0 * exp(η.γ), βγ => θ.βγ,
        α0 => θ.α0 * exp(η.α), βα => θ.βα, δ => θ.δ * exp(η.δ), n => θ.n, AUC => c.AUC],
    init = (θ, η, c) -> [TUMOUR.P => θ.P0, TUMOUR.K => θ.K0 * exp(η.K0)],
    observed = (P = TUMOUR.P,),
    error = (P = CombinedError(3.0, 0.14),),
    tspan = (0.0, 75.0), solver = Tsit5(), solve_kwargs = TOL)

days = collect(0.0:1.0:75.0)
nomouse = Regimen(Dose[])                   # empty regimen
mice = virtual_population(500; regimen = nomouse, obstimes = days, seed = 2015,
    covariate_sampler = rng -> (AUC = MOUSE_AUC.dose100,))
res = simulate(model, mice; seed = 2015)
println(res)

# ---- check: zero-random-effect prediction equals the atlas typical curve ----------------
typical = zero_re_prediction(model, mice[1])
prob, _ = problem_for("fig4_left_P"; fixes = FIXES)          # atlas "with fixes" mouse panel: P0 150, AUC 1140.8
atlas = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10, saveat = days)
dev = maximum(abs.([o.ipred for o in typical.obs] .- atlas[TUMOUR.P]))
println("max |zero-RE prediction - atlas typical curve| = ", dev, " mm^3 (curve spans ",
    round(minimum(atlas[TUMOUR.P]); digits = 1), " to ", round(maximum(atlas[TUMOUR.P]); digits = 1), ")")
@assert dev < 1e-6

# ---- percentile bands -------------------------------------------------------------------
function band(field)
    q = [quantile([getfield(o, field) for o in res.obs if o.time == d], [0.05, 0.5, 0.95]) for d in days]
    return (lo = [x[1] for x in q], med = [x[2] for x in q], hi = [x[3] for x in q])
end
ip, dv = band(:ipred), band(:dv)

M = CairoMakie
fig = M.Figure(size = (760, 480))
ax = M.Axis(fig[1, 1]; xlabel = "Time (days)", ylabel = "Tumour volume (mm³)",
    title = "Ouerdani 2015 mice, 100 mg/kg: 500 virtual mice")
M.band!(ax, days, dv.lo, dv.hi; color = (:steelblue, 0.2), label = "5-95 % observed (with residual error)")
M.band!(ax, days, ip.lo, ip.hi; color = (:steelblue, 0.45), label = "5-95 % individual prediction")
M.lines!(ax, days, ip.med; color = :steelblue, linewidth = 2.5, label = "median")
M.lines!(ax, days, [o.ipred for o in typical.obs]; color = :black, linestyle = :dash, linewidth = 1.5, label = "zero-random-effect (atlas typical)")
M.axislegend(ax; position = :rt)
out = joinpath(@__DIR__, "population_ouerdani.png")
M.save(out, fig)
println("saved ", out)
