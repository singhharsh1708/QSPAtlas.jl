# Population simulation performance.
#
#   julia --project=benchmarks/solvers --threads=auto benchmarks/population/run.jl
#
# Times `simulate` for cohorts of 100, 1,000 and 10,000 virtual subjects on two models, and
# one simulation-based VPC. Compilation is timed separately (first call); the reported times
# are the median of 3 runs after that. Results: results/summary.md and results/summary.toml.
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using Statistics: median
using Printf, TOML

const OUT = joinpath(@__DIR__, "results")
mkpath(OUT)

# ---- model 1: Ouerdani 2015 tumour growth / angiogenesis, mice, 500 -> N virtual mice -------
Base.include(@__MODULE__, joinpath(models_dir(), "ouerdani_2015", "model.jl"))
Base.include(@__MODULE__, joinpath(models_dir(), "ouerdani_2015", "protocol.jl"))

ηn = [:λ, :K0, :γ, :α, :δ]
ω = (λ = 0.53, K0 = 0.36, γ = 0.19, α = 0.24, δ = 0.42)                  # Table 1 IIV
ouerdani = PopulationModel(TUMOUR;
    routes = Route[], θ = (PRECLINICAL..., P0 = FIXED_P0.mouse), ηnames = ηn,
    Ω = [ω[a]^2 * (i == j) for (i, a) in enumerate(ηn), (j, _) in enumerate(ηn)],
    params = (θ, η, c) -> [λ => θ.λ * exp(η.λ), b => θ.b, γ0 => θ.γ0 * exp(η.γ), βγ => θ.βγ,
        α0 => θ.α0 * exp(η.α), βα => θ.βα, δ => θ.δ * exp(η.δ), n => θ.n, AUC => c.AUC],
    init = (θ, η, c) -> [TUMOUR.P => θ.P0, TUMOUR.K => θ.K0 * exp(η.K0)],
    observed = (P = TUMOUR.P,), error = (P = CombinedError(3.0, 0.14),),
    tspan = (0.0, 75.0), solver = Tsit5())
ouerdani_cohort(n) = virtual_population(n; regimen = Regimen(Dose[]), obstimes = collect(0.0:1.0:75.0),
    seed = 1, covariate_sampler = rng -> (AUC = MOUSE_AUC.dose100,))

# ---- model 2: two-compartment PK, 1 h infusion every 12 h for 14 days, weight covariate ------
@parameters CL = 5.0 V1 = 20.0 Q = 10.0 V2 = 40.0 Rin = 0.0
@variables A1(t) = 0.0 A2(t) = 0.0
@mtkcompile pk2 = System([D(A1) ~ -CL / V1 * A1 - Q / V1 * A1 + Q / V2 * A2 + Rin,
                          D(A2) ~ Q / V1 * A1 - Q / V2 * A2], t)
twocmt = PopulationModel(pk2;
    routes = [Route(:central, pk2.A1; rate = pk2.Rin)],
    θ = (CL = 5.0, V1 = 20.0, Q = 10.0, V2 = 40.0), ηnames = [:CL, :V1], Ω = [0.09 0.02; 0.02 0.04],
    params = (θ, η, c) -> [pk2.CL => θ.CL * exp(η.CL) * (c.WT / 70)^0.75, pk2.V1 => θ.V1 * exp(η.V1) * (c.WT / 70),
        pk2.Q => θ.Q, pk2.V2 => θ.V2],
    observed = (conc = pk2.A1 / pk2.V1,), error = (conc = ProportionalError(0.15),),
    tspan = (0.0, 336.0), solver = Tsit5())
twocmt_cohort(n) = virtual_population(n; regimen = Regimen(100.0; cmt = :central, duration = 1.0, ii = 12.0, addl = 27),
    obstimes = collect(0.0:2.0:336.0), seed = 1, covariate_sampler = rng -> (WT = 50 + 50rand(rng),))

# ---- timing ----------------------------------------------------------------------------------
function time_cohort(model, cohort, n)
    subs = cohort(n)
    times = [(@elapsed simulate(model, subs; seed = 1)) for _ in 1:3]
    return median(times)
end

rows = []
for (name, model, cohort) in (("Ouerdani 2015 tumour model (2 states, no dosing)", ouerdani, ouerdani_cohort),
                              ("Two-compartment PK, 28 infusions, weight covariate", twocmt, twocmt_cohort))
    tc = @elapsed simulate(model, cohort(10); seed = 1)              # includes compilation
    for n in (100, 1_000, 10_000)
        t = time_cohort(model, cohort, n)
        push!(rows, (; name, n, seconds = t, per_subject_ms = 1e3 * t / n, compile = tc))
        @info "timed" name n t
    end
end

vsubs = twocmt_cohort(100)
vpc(twocmt, vsubs[1:10]; endpoint = :conc, bins = collect(0.0:24.0:336.0), replicates = 10)   # compile
tvpc = @elapsed vpc(twocmt, vsubs; endpoint = :conc, q = (0.05, 0.5, 0.95), bins = collect(0.0:24.0:336.0), replicates = 200)

# ---- report ----------------------------------------------------------------------------------
cpu = Sys.cpu_info()[1].model
io = IOBuffer()
println(io, "Machine: ", cpu, ", ", Threads.nthreads(), " Julia threads, Julia ", VERSION, ".\n")
println(io, "| model | subjects | wall time | per subject |")
println(io, "|---|---|---|---|")
for r in rows
    println(io, "| ", r.name, " | ", r.n, " | ", @sprintf("%.2f s", r.seconds), " | ", @sprintf("%.2f ms", r.per_subject_ms), " |")
end
println(io, "\nVPC of the two-compartment model, 100 subjects x 200 replicates (20,000 simulated subjects), ",
    "3 quantiles in 14 bins: ", @sprintf("%.1f s", tvpc), ".\n")
println(io, "First-call compile times: ", join([@sprintf("%s %.0f s", r.name, r.compile) for r in rows if r.n == 100], "; "), ".")
write(joinpath(OUT, "summary.md"), String(take!(io)))
open(f -> TOML.print(f, Dict("cpu" => cpu, "threads" => Threads.nthreads(), "julia" => string(VERSION), "vpc_seconds" => tvpc,
    "rows" => [Dict(string(k) => v for (k, v) in pairs(r)) for r in rows])), joinpath(OUT, "summary.toml"), "w")
print(read(joinpath(OUT, "summary.md"), String))
