# One-at-a-time sensitivity of the Figure 3/4 comparison to inputs that the paper does not pin down.
# INFORMATIONAL ONLY: nothing here is a fix, grade input or fitted value. The readout fix (unbound-only
# periphery/tumour concentrations) is applied in every variant so that the other inputs can be seen.
#   julia --project=. models/stroh_2019/sensitivity.jl > sensitivity.log
using QSPAtlas, Printf
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
Base.include(Main, joinpath(@__DIR__, "model.jl"))
Base.include(Main, joinpath(@__DIR__, "protocol.jl"))

refs = sort(collect(values(load_references(joinpath(@__DIR__, "reference")))); by = r -> r.id)
variants = Dict{String, Function}(
    "base (kcleave 3e-7, tumour kcleave x10, kopen fixed)" => (sp, p) -> p,
    "kcleave x0.33 (1e-7)" => (sp, p) -> (for k in (:kcleave_central, :kcleave_peripheral, :kcleave_tumor); p[k] *= 1 / 3; end; p),
    "kcleave x0 (0)" => (sp, p) -> (for k in (:kcleave_central, :kcleave_peripheral, :kcleave_tumor); p[k] = 0.0; end; p),
    "tumour kcleave x1 (no tumour enhancement)" => (sp, p) -> (p[:kcleave_tumor] = p[:kcleave_central]; p),
    "ksynR tumour 2.8221e-5 (supplement parameter 28) instead of 2.7e-4" => (sp, p) -> (p[:ksyn_tumor] = 2.8221e-5; p),
    "kclose = 0.232 fixed, kopen = kclose/Kmask" => (sp, p) -> (p[:kopen] = 0.232 / (p[:kclose] / p[:kopen]); p[:kclose] = 0.232; p),
    "k13, k31 as supplement (2.5e-9, 4.4e-6)" => (sp, p) -> (p[:k13] = 2.5e-9; p[:k31] = 4.4e-6; p),
)
function sim(sp, out, mod)
    sys = get_system()
    pars, dosed, bw = molecule_setup(sp; fixes = FIXES)
    pars = mod(sp, copy(pars))
    prob = ODEProblem(sys, problem_pairs(sys, pars), (0.0, sp.tend))
    cb = bolus_callback(sys, getproperty(sys, dosed), [0.0], dose_nM(sp.dose, bw, pars[:V1]))
    ts = collect(range(0.0, sp.tend; length = 1500))
    sol = solve(prob, Rodas5P(); callback = cb, saveat = ts, abstol = 1e-10, reltol = 1e-8)
    return ts, sol(ts; idxs = getproperty(sys, out)).u
end
for (name, mod) in sort(collect(variants); by = first)
    res = Dict("fig3" => Float64[], "fig4" => Float64[]); npass = Dict("fig3" => 0, "fig4" => 0)
    cache = Dict{Any, Any}()
    for ref in refs
        sp, out = spec_for(ref.id)
        out = Symbol(replace(String(out), "Ctot_" => "Cfree_"))
        ts, y = get!(cache, (sp, out)) do; sim(sp, out, mod); end
        c = compare_curve(ref, ts, y; tolerance = TOLERANCE)
        g = startswith(ref.id, "fig3") ? "fig3" : "fig4"
        push!(res[g], c.rms); npass[g] += c.passed
    end
    @printf("%-75s fig3 pass %2d/%d mean rms %.3f | fig4 pass %2d/%d mean rms %.3f\n", name, npass["fig3"], length(res["fig3"]), sum(res["fig3"]) / length(res["fig3"]),
        npass["fig4"], length(res["fig4"]), sum(res["fig4"]) / length(res["fig4"]))
end
