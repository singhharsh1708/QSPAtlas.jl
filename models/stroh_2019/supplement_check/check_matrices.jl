# Structural check of model.jl against the supplement matrices (PSP4-8-676-s003.pdf, transcribed to A1.csv / A2.csv).
# Builds A1 and A2 from the reaction list of model.jl with the parameter values of the supplement's printed
# parameter list (k.csv) and unit volume ratios, and compares entry by entry.
#   julia --project=. models/stroh_2019/supplement_check/check_matrices.jl
using DelimitedFiles
Base.include(Main, joinpath(@__DIR__, "..", "model.jl"))

# supplement state order (N = 55)
idx = Dict{Symbol, Int}()
let n = 0
    for comp in COMPS
        idx[rname(comp)] = (n += 1)
        for s in SPEC; idx[sname(s, comp)] = (n += 1); end
        if comp === :central
            idx[:pb_el_central] = (n += 1)
        else
            idx[Symbol("pb_tmdd_m_", comp)] = (n += 1); idx[Symbol("pb_tmdd_o_", comp)] = (n += 1)
            idx[Symbol("pb_el_", comp)] = (n += 1)
        end
    end
end
N = length(idx)
@assert N == 55 && idx[:R_peripheral] == 18 && idx[:R_tumor] == 37 && idx[:pb_el_peripheral] == 36

pv = Dict(:kopen => 1.16e-2, :kclose => 0.232, :kcleave_central => 1e-6, :kcleave_peripheral => 1e-6, :kcleave_tumor => 1e-5,
    :kel => 6.3e-7, :kon1 => 1e-3, :kon2 => 1e-3, :koff1 => 2e-3, :kendo => 1e-4, :k12 => 4.8e-6, :k21 => 4.4e-6, :k31 => 4.4e-6, :k13 => 2.5e-9)
A1 = Dict{Tuple{Int, Int}, Float64}(); A2 = Dict{Tuple{Int, Int}, Float64}()
add!(D, k, v) = (D[k] = get(D, k, 0.0) + v)
for rx in reactions()
    r = pv[rx.rate] * rx.mult
    if length(rx.from) == 1
        j = idx[rx.from[1]]
        add!(A1, (j, j), -r)
        for s in rx.to; add!(A1, (idx[s], j), r); end
    else
        i, j = idx[rx.from[1]], idx[rx.from[2]]
        col = (i - 1) * N + j
        add!(A2, (i, col), -r); add!(A2, (j, col), -r)
        for s in rx.to; add!(A2, (idx[s], col), r); end
    end
end
function load(f)
    d = readdlm(joinpath(@__DIR__, f), ','; skipstart = 1)
    D = Dict{Tuple{Int, Int}, Float64}()
    for k in 1:size(d, 1); add!(D, (Int(d[k, 1]), Int(d[k, 2])), Float64(d[k, 3])); end
    return D
end
function compare(mine, sup, name)
    bad = 0
    for k in union(keys(mine), keys(sup))
        a, b = get(mine, k, 0.0), get(sup, k, 0.0)
        if abs(a - b) > 1e-9 * max(abs(a), abs(b), 1e-12) + 1e-18
            bad += 1
            bad <= 25 && println("  $name mismatch at $k: model $a, supplement $b")
        end
    end
    println(name, ": ", length(sup), " supplement entries, ", length(mine), " model entries, ", bad, " mismatches")
    return bad
end
bad = compare(A1, load("A1.csv"), "A1") + compare(A2, load("A2.csv"), "A2")
println(bad == 0 ? "ALL ENTRIES MATCH" : "MISMATCHES: $bad")
