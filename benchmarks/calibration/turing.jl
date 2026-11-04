# Turing.jl (NUTS) recalibration of hutchinson_2016 against Figure 1d (group means).
# Run from QSPAtlas.jl/:  julia -t 4 --startup-file=no --project=benchmarks/calibration benchmarks/calibration/turing.jl
# Env: NSAMP (default 1000), NCHAINS (4), NADAPT (500), PRIOR = published | generic
using Dates
println("start ", now(), " threads ", Threads.nthreads())
include(joinpath(@__DIR__, "common.jl"))
using Turing, Distributions, Random, MCMCChains, StatsBase, LinearAlgebra, TOML, Statistics, ForwardDiff
using DynamicPPL: InitFromParams

NSAMP = parse(Int, get(ENV, "NSAMP", "1000"))
NCH = parse(Int, get(ENV, "NCHAINS", "4"))
NADAPT = parse(Int, get(ENV, "NADAPT", "500"))
PRIOR = get(ENV, "PRIOR", "published")
tag = PRIOR == "published" ? "" : "_" * PRIOR

# prior medians: published Table 1 values (sd 1 on the log scale), or a generic set not taken from the paper (sd 1.5)
pm = PRIOR == "published" ?
    (alpha_T = PUB.alpha_T, alpha_V = PUB.alpha_V, K = PUB.K, Nmax = PUB.Nmax, delta_V_bev = PUB.delta_V_bev, delta_V_van = PUB.delta_V_van, sigma = 0.2) :
    (alpha_T = 0.1, alpha_V = 0.1, K = 1.0, Nmax = 3.0, delta_V_bev = 0.05, delta_V_van = 0.05, sigma = 0.2)
psd = PRIOR == "published" ? 1.0 : 1.5

df = load_data()
ts = Dict(g => Float64.(df.time[df.group .== g]) for g in GROUPS)
ylog = Dict(g => log.(df.volume[df.group .== g]) for g in GROUPS)
prob_sim = make_prob()

@model function hutch(prob, ts, ylog, pm, psd)
    alpha_T ~ LogNormal(log(pm.alpha_T), psd)
    alpha_V ~ LogNormal(log(pm.alpha_V), psd)
    K ~ LogNormal(log(pm.K), psd)
    Nmax ~ truncated(LogNormal(log(pm.Nmax), psd); lower = 1.0)
    delta_V_bev ~ LogNormal(log(pm.delta_V_bev), psd)
    delta_V_van ~ LogNormal(log(pm.delta_V_van), psd)
    sigma ~ LogNormal(log(pm.sigma), psd)
    θ = (; alpha_T, alpha_V, K, Nmax, delta_V_bev, delta_V_van)
    for g in GROUPS
        yp = simulate_group(prob, g, θ, ts[g])
        # log-normal residual: log y ~ Normal(log T, sigma)   (the -log y Jacobian is constant and omitted)
        for i in eachindex(yp)
            Turing.@addlogprob! logpdf(Normal(log(yp[i]), sigma), ylog[g][i])
        end
    end
end
model = hutch(prob_sim, ts, ylog, pm, psd)

rng = MersenneTwister(20261104)
inits = [InitFromParams((; (k => (k == :Nmax ? max(1.05, pm[k] * exp(0.4 * randn(rng))) : pm[k] * exp(0.4 * randn(rng))) for k in keys(pm))...)) for _ in 1:NCH]

t_s = @elapsed chn = sample(rng, model, NUTS(NADAPT, 0.8), MCMCThreads(), NSAMP, NCH; initial_params = inits, progress = false, chain_type = MCMCChains.Chains)
println("sampling done in ", round(t_s; digits = 1), " s")

ss = DataFrame(summarystats(chn))
qs = DataFrame(quantile(chn; q = [0.05, 0.5, 0.95]))
println(ss); println(qs)
# divergences and tree depth
nm = chn[:, [:numerical_error, :tree_depth, :n_steps, :acceptance_rate, :step_size], :]
println("divergent transitions: ", sum(Array(chn[:, :numerical_error, :])), " of ", length(chn))
println("mean tree depth ", mean(Array(chn[:, :tree_depth, :])), "  mean step size ", mean(Array(chn[:, :step_size, :])), "  mean accept ", mean(Array(chn[:, :acceptance_rate, :])))

# draws
pn = collect(keys(pm))
draws = DataFrame(; (k => vec(Array(chn[:, k, :])) for k in pn)...)
draws.chain = repeat(1:NCH; inner = NSAMP)
CSV.write(joinpath(@__DIR__, "results", "turing_draws$tag.csv"), draws)
res = DataFrame(parameter = String[], median = Float64[], q05 = Float64[], q95 = Float64[], mean = Float64[], rhat = Float64[], ess_bulk = Float64[], ess_tail = Float64[])
for k in pn
    x = draws[!, k]
    r = ss[ss.parameters .== k, :]
    push!(res, (String(k), median(x), quantile(x, 0.05), quantile(x, 0.95), mean(x), r.rhat[1], r.ess_bulk[1], r.ess_tail[1]))
end
CSV.write(joinpath(@__DIR__, "results", "turing_summary$tag.csv"), res)
println(res)

# pairwise correlations on the log scale
L = log.(Matrix(draws[:, pn]))
C = cor(L)
cdf = DataFrame(C, String.(pn)); insertcols!(cdf, 1, :parameter => String.(pn))
CSV.write(joinpath(@__DIR__, "results", "turing_logcorr$tag.csv"), cdf)
show(stdout, "text/plain", round.(C; digits = 2)); println()

# posterior predictive: typical-individual curve (no noise) and with log-normal noise on the data days
rng2 = MersenneTwister(7)
idx = rand(rng2, 1:nrow(draws), 400)
pp = DataFrame(group = String[], time = Float64[], q05 = Float64[], q50 = Float64[], q95 = Float64[])
curves = Dict(g => zeros(length(idx), length(TGRID)) for g in GROUPS)
for (j, i) in enumerate(idx)
    θ = (alpha_T = draws.alpha_T[i], alpha_V = draws.alpha_V[i], K = draws.K[i], Nmax = draws.Nmax[i], delta_V_bev = draws.delta_V_bev[i], delta_V_van = draws.delta_V_van[i])
    for g in GROUPS
        curves[g][j, :] = simulate_group(prob_sim, g, θ, TGRID)
    end
end
for g in GROUPS, (k, tt) in enumerate(TGRID)
    v = curves[g][:, k]
    push!(pp, (g, tt, quantile(v, 0.05), quantile(v, 0.5), quantile(v, 0.95)))
end
CSV.write(joinpath(@__DIR__, "results", "turing_predictive$tag.csv"), pp)

# lp at posterior draws vs published: maximum a posteriori / best draw
println("done ", now())
