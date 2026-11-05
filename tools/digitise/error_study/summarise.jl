# Summaries of the error study, the tolerance rule, and its effect on the committed scores.
#   julia --project=tools/digitise tools/digitise/error_study/summarise.jl
# Reads results/synthetic_runs.csv (synth.jl), results/real_figures.csv (real_figures/compare_real.jl)
# and the per-curve rms/p95 in models/*/score.toml (not rescored). Writes results/synthetic_summary.csv,
# results/tolerance_changes.csv and prints the tables used in README.md.
using DelimitedFiles, Statistics, Printf, TOML

const RES = joinpath(@__DIR__, "results")
const MODELS = joinpath(@__DIR__, "..", "..", "..", "models")

function readcsv(path)
    m, h = readdlm(path, ','; header = true)
    return Dict(String(h[k]) => m[:, k] for k in eachindex(h))
end
q(v, p) = isempty(v) ? NaN : quantile(v, p)
f4(x) = @sprintf("%.4f", x)

function synthetic_summary()
    d = readcsv(joinpath(RES, "synthetic_runs.csv"))
    panels = unique(d["panel"])
    rows = String[]
    println("\n| panel type | tracing rule | curves | not traced | RMS median / 95th / 99th | P95 median / 95th / 99th | RMS > 0.015 | fails 0.015/0.03 |")
    println("|---|---|---|---|---|---|---|---|")
    allr = Float64[]; allp = Float64[]
    for p in panels
        k = d["panel"] .== p
        ok = k .& (d["traced"] .== true)
        r = Float64.(d["rms"][ok]); s = Float64.(d["p95"][ok])
        r0 = Float64.(d["rms_trace_only"][ok]); s0 = Float64.(d["p95_trace_only"][ok])
        append!(allr, r); append!(allp, s)
        nf = count(k) - count(ok)
        fail = count((r .> 0.015) .| (s .> 0.03))
        push!(rows, join((p, d["rule"][findfirst(k)], count(k), nf, f4(median(r)), f4(q(r, 0.95)), f4(q(r, 0.99)), f4(maximum(r)),
                          f4(median(s)), f4(q(s, 0.95)), f4(q(s, 0.99)), f4(maximum(s)), f4(q(r0, 0.99)), f4(q(s0, 0.99)),
                          f4(mean(r .> 0.015)), f4(fail / length(r)), f4(median(Float64.(d["calerr_px"][k])))), ","))
        @printf("| %s | %s | %d | %d | %.4f / %.4f / %.4f | %.4f / %.4f / %.4f | %.1f%% | %.1f%% |\n", p, d["rule"][findfirst(k)], count(k), nf,
                median(r), q(r, 0.95), q(r, 0.99), median(s), q(s, 0.95), q(s, 0.99), 100mean(r .> 0.015), 100fail / length(r))
    end
    @printf("| all | | %d | | %.4f / %.4f / %.4f | %.4f / %.4f / %.4f | %.1f%% | %.1f%% |\n", length(allr), median(allr), q(allr, 0.95), q(allr, 0.99),
            median(allp), q(allp, 0.95), q(allp, 0.99), 100mean(allr .> 0.015), 100mean((allr .> 0.015) .| (allp .> 0.03)))
    open(joinpath(RES, "synthetic_summary.csv"), "w") do io
        println(io, "panel,rule,curves,not_traced,rms_median,rms_q95,rms_q99,rms_max,p95_median,p95_q95,p95_q99,p95_max,rms_q99_trace_only,p95_q99_trace_only,frac_rms_over_0.015,frac_fail_current_tol,calibration_err_px_median")
        foreach(l -> println(io, l), rows)
        println(io, join(("all", "", length(allr), "", f4(median(allr)), f4(q(allr, 0.95)), f4(q(allr, 0.99)), f4(maximum(allr)), f4(median(allp)),
                          f4(q(allp, 0.95)), f4(q(allp, 0.99)), f4(maximum(allp)), "", "", f4(mean(allr .> 0.015)), f4(mean((allr .> 0.015) .| (allp .> 0.03))), ""), ","))
    end
    return d
end

function tolerance_changes(tr, tp)
    out = String[]
    changed = String[]
    counts = Dict{Tuple{String, String}, NTuple{3, Int}}()
    for m in sort(filter(x -> isfile(joinpath(MODELS, x, "score.toml")), readdir(MODELS)))
        s = TOML.parsefile(joinpath(MODELS, m, "score.toml"))
        for run in ("published", "fixed")
            haskey(s, run) || continue
            n = 0; before = 0; after = 0
            for c in get(s[run], "curves", [])
                rms, p95 = c["rms"], c["p95"]
                old = c["passed"]
                new = rms <= tr && p95 <= tp
                cur = rms <= c["tol_rms"] && p95 <= c["tol_p95"]
                cur == old || @warn "stored pass flag differs from its own tolerance" m c["id"]
                n += 1; before += old; after += new
                push!(out, join((m, run, c["id"], rms, p95, old, new), ","))
                old != new && push!(changed, @sprintf("%s %s %s: rms %.4f p95 %.4f, %s -> %s", m, run, c["id"], rms, p95, old ? "pass" : "fail", new ? "pass" : "fail"))
            end
            counts[(m, run)] = (n, before, after)
        end
    end
    open(joinpath(RES, "tolerance_changes.csv"), "w") do io
        println(io, "model,run,curve_id,rms,p95,pass_current_0.015_0.03,pass_new_$(tr)_$(tp)")
        foreach(l -> println(io, l), out)
    end
    @printf("\nTolerance rms <= %.4f and p95 <= %.4f against the committed scores\n", tr, tp)
    println("| model | run | curves | pass now | pass at new tolerance |")
    println("|---|---|---|---|---|")
    for k in sort(collect(keys(counts)))
        n, b, a = counts[k]; println("| $(k[1]) | $(k[2]) | $n | $b | $a |")
    end
    println("\nchanged:"); foreach(println, changed)
end

function real_summary()
    f = joinpath(RES, "real_figures.csv")
    isfile(f) || return
    d = readcsv(f)
    println("\nReal figures: atlas reference vs thesis curve")
    for m in unique(d["model"])
        k = d["model"] .== m
        @printf("%s: %d curves, rms median %.4f max %.4f, p95 median %.4f max %.4f\n", m, count(k), median(Float64.(d["rms"][k])), maximum(Float64.(d["rms"][k])),
                median(Float64.(d["p95"][k])), maximum(Float64.(d["p95"][k])))
    end
end

# Overlay check: the atlas accepts a trace only after its points are drawn over the figure. A point more
# than ACCEPT of the axis box off the line (4 to 15 px on the atlas panels) is plainly visible there, and
# such traces are fixed (exclusions, outlier filters) before use. Synthetic traces with a larger maximum
# are counted as caught by the check and left out of the tolerance distribution.
const ACCEPT = 0.05
const SAFETY = 2.0

function accepted_summary()
    d = readcsv(joinpath(RES, "synthetic_runs.csv"))
    ok = (d["traced"] .== true) .& (Float64.(d["max"]) .<= ACCEPT)
    out = Dict{String, NTuple{2, Float64}}()
    rows = String[]
    println("\nAccepted traces (max point distance <= $ACCEPT)")
    println("| panel type | accepted | caught by the check | RMS median / 95th / 99th | P95 median / 95th / 99th |")
    println("|---|---|---|---|---|")
    for p in unique(d["panel"])
        k = (d["panel"] .== p) .& (d["traced"] .== true); a = k .& ok
        r = Float64.(d["rms"][a]); s = Float64.(d["p95"][a])
        out[p] = (q(r, 0.99), q(s, 0.99))
        @printf("| %s | %d | %.1f%% | %.4f / %.4f / %.4f | %.4f / %.4f / %.4f |\n", p, count(a), 100(1 - count(a) / count(k)),
                median(r), q(r, 0.95), q(r, 0.99), median(s), q(s, 0.95), q(s, 0.99))
        push!(rows, join((p, count(a), f4(1 - count(a) / count(k)), f4(median(r)), f4(q(r, 0.95)), f4(q(r, 0.99)), f4(maximum(r)),
                          f4(median(s)), f4(q(s, 0.95)), f4(q(s, 0.99)), f4(maximum(s))), ","))
    end
    r = Float64.(d["rms"][ok]); s = Float64.(d["p95"][ok])
    @printf("| all | %d | %.1f%% | %.4f / %.4f / %.4f | %.4f / %.4f / %.4f |\n", count(ok), 100(1 - count(ok) / count(d["traced"] .== true)),
            median(r), q(r, 0.95), q(r, 0.99), median(s), q(s, 0.95), q(s, 0.99))
    open(joinpath(RES, "synthetic_summary_accepted.csv"), "w") do io
        println(io, "panel,accepted,frac_caught_by_overlay_check,rms_median,rms_q95,rms_q99,rms_max,p95_median,p95_q95,p95_q99,p95_max")
        foreach(l -> println(io, l), rows)
    end
    return out
end

# panel type of a committed curve (for the per-type alternative)
function panel_type(m, id)
    m == "hutchinson_2016" && return startswith(id, "fig6") ? "hut_fig6" : (startswith(id, "fig2_a") || startswith(id, "fig2_b")) ? "hut_fig2_single" : "hut_fig2_five"
    m == "ouerdani_2015" && return "ouerdani_fig4"
    m == "stroh_2019" && return startswith(id, "fig3") ? "stroh_fig3" : (startswith(id, "fig4d") || startswith(id, "fig4e")) ? "stroh_fig4_lin" : "stroh_fig4_log"
    m == "fornari_2019" && return startswith(id, "fig3") ? "fornari_fig3" : "fornari_fig4"
    return ""
end

function flips(name, tolfor)
    n = 0; changed = String[]
    open(joinpath(RES, "tolerance_changes_$(name).csv"), "w") do io
        println(io, "model,run,curve_id,rms,p95,tol_rms,tol_p95,pass_current,pass_new")
        for m in sort(filter(x -> isfile(joinpath(MODELS, x, "score.toml")), readdir(MODELS)))
            sc = TOML.parsefile(joinpath(MODELS, m, "score.toml"))
            for run in ("published", "fixed"), c in get(get(sc, run, Dict()), "curves", [])
                tr, tp = tolfor(m, c["id"])
                new = c["rms"] <= tr && c["p95"] <= tp
                println(io, join((m, run, c["id"], c["rms"], c["p95"], tr, tp, c["passed"], new), ","))
                n += 1
                new != c["passed"] && push!(changed, @sprintf("%s %s %s (rms %.4f, p95 %.4f): %s", m, run, c["id"], c["rms"], c["p95"], new ? "fail -> pass" : "pass -> fail"))
            end
        end
    end
    println("\n$name: $(length(changed)) of $n committed curve results change")
    foreach(x -> println("  ", x), changed)
end

synthetic_summary()
real_summary()
q99 = accepted_summary()
worst_r = maximum(first.(values(q99))); worst_p = maximum(last.(values(q99)))
tr = round(SAFETY * worst_r; sigdigits = 2); tp = round(SAFETY * worst_p; sigdigits = 2)
@printf("\nRule: tolerance = %.0f x the largest per-panel-type 99th percentile of accepted digitising error\n", SAFETY)
@printf("  RMS: %.0f x %.4f = %.4f -> %.3f;  P95: %.0f x %.4f = %.4f -> %.3f\n", SAFETY, worst_r, SAFETY * worst_r, tr, SAFETY, worst_p, SAFETY * worst_p, tp)
open(joinpath(RES, "recommended_tolerance.toml"), "w") do io
    TOML.print(io, Dict("rule" => "tolerance = $(SAFETY) x largest per-panel-type 99th percentile of accepted synthetic digitising error, 2 significant figures",
                        "accept_max_point_distance" => ACCEPT, "safety_factor" => SAFETY, "worst_q99_rms" => worst_r, "worst_q99_p95" => worst_p,
                        "rms" => tr, "p95" => tp))
end
flips("recommended", (m, id) -> (tr, tp))
flips("per_panel_type", (m, id) -> (t = panel_type(m, id); t == "" ? (tr, tp) : (round(SAFETY * q99[t][1]; sigdigits = 2), round(SAFETY * q99[t][2]; sigdigits = 2))))
