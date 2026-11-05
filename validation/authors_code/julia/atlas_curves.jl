# Atlas curves for the authors'-code comparison. Writes out/atlas_<model>.csv (case, time, var, value).
# Each model's model.jl and protocol.jl are loaded into their own module; the atlas functions are called
# with the atlas's documented fixes (FIXES), so the curves are the ones the atlas grades.
# Run from QSPAtlas.jl: julia --startup-file=no --project=. validation/authors_code/julia/atlas_curves.jl
using QSPAtlas, ModelingToolkit, OrdinaryDiffEq, DelimitedFiles

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUT = joinpath(@__DIR__, "..", "out")
mkpath(OUT)

function load_model(id)
    m = Module(Symbol(id))
    Core.eval(m, :(using QSPAtlas, ModelingToolkit, OrdinaryDiffEq))
    Base.include(m, joinpath(ROOT, "models", id, "model.jl"))
    Base.include(m, joinpath(ROOT, "models", id, "protocol.jl"))
    return m
end

function write_long(path, rows)
    open(path, "w") do io
        println(io, "case,time,var,value")
        for (c, t, v, y) in rows
            println(io, c, ",", repr(Float64(t)), ",", v, ",", repr(Float64(y)))
        end
    end
    println("wrote ", path, " (", length(rows), " rows)")
end

# ---------------------------------------------------------------- Ouerdani 2015 (Figure 4 typical)
let M = load_model("ouerdani_2015")
    rows = Tuple{String, Float64, String, Float64}[]
    sys = M.TUMOUR
    # Figure 4 cases with the atlas fixes; extra mouse doses (10, 30 mg/kg AUC) as structure checks
    cases = [("mouse_fig4", "fig4_left_P", nothing, 0:0.25:75),
             ("mouse_auc220", "fig4_left_P", M.MOUSE_AUC.dose10, 0:0.25:75),
             ("mouse_auc657", "fig4_left_P", M.MOUSE_AUC.dose30, 0:0.25:75),
             ("patient_fig4", "fig4_right_P", nothing, 0:5.0:2500),
             ("mouse_control_auc0", "fig4_left_P", 0.0, 0:0.25:75)]
    for (name, key, auc, grid) in cases
        prob, _ = M.problem_for(key; fixes = M.FIXES, AUC = auc)
        prob = remake(prob; tspan = (0.0, last(grid)))
        sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10, saveat = collect(grid))
        for (i, tt) in enumerate(sol.t)
            push!(rows, (name, tt, "P", sol[sys.P][i]))
            push!(rows, (name, tt, "K", sol[sys.K][i]))
        end
    end
    write_long(joinpath(OUT, "atlas_ouerdani_2015.csv"), rows)
end

# ---------------------------------------------------------------- Hutchinson 2016 (Figure 2, Figure 6a)
let M = load_model("hutchinson_2016")
    rows = Tuple{String, Float64, String, Float64}[]
    ids = ["fig2_a_T", "fig2_c_T_d04", "fig2_c_T_d06", "fig2_c_T_d08", "fig2_c_T_d10", "fig2_c_T_d12",
           "fig2_e_T_N2", "fig2_e_T_N3", "fig2_e_T_N6", "fig2_e_T_N11", "fig2_e_T_N21", "fig6_a_T"]
    for id in ids
        p, _ = M.scenario(id)
        s = M.simulate_scenario(p; fixes = M.FIXES)
        case = replace(id, "_T" => "")
        for i in eachindex(s.t)
            push!(rows, (case, s.t[i], "T", s.T[i]))
            push!(rows, (case, s.t[i], "V", s.V[i]))
            push!(rows, (case, s.t[i], "N", s.N[i]))
        end
    end
    write_long(joinpath(OUT, "atlas_hutchinson_2016.csv"), rows)
end

# ---------------------------------------------------------------- Stein 2019 (typical value curve)
let M = load_model("stein_2019")
    rows = Tuple{String, Float64, String, Float64}[]
    th = M.TYPICAL
    grid = sort(unique(vcat(collect(0:0.1:30), collect(31:1.0:380))))
    # (case, toci day, steroid day); Inf = not given
    cases = [("typical", Inf, Inf), ("toci4_ster6", 4.0, 6.0), ("ster2_toci6", 6.0, 2.0), ("toci12", 12.0, Inf)]
    for (name, toci, ster) in cases
        Tm = th.Tmax
        if toci <= ster
            t1, f1, t2, f2 = min(toci, Tm), th.Ftoci, min(ster, Tm), th.Fster
        else
            t1, f1, t2, f2 = min(ster, Tm), th.Fster, min(toci, Tm), th.Ftoci
        end
        sys = M.CART
        prob = ODEProblem(sys, [sys.E => th.Cmax / th.foldx, sys.M => 0.0, sys.foldx => th.foldx, sys.Tmax => Tm,
            sys.α => th.α, sys.FB => th.FB, sys.β => th.β, sys.F1 => f1, sys.F2 => f2, sys.T1 => t1, sys.T2 => t2],
            (0.0, 380.0))
        sol = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12, saveat = grid, tstops = [t1, t2, Tm])
        C = sol[sys.C]
        for i in eachindex(sol.t)
            push!(rows, (name, sol.t[i], "C", C[i]))
            # closed form of the atlas (Eq. 1/2), for reference
            push!(rows, (name, sol.t[i], "C_closed", M.closed_form(sol.t[i], (foldx = th.foldx, Tmax = Tm, Cmax = th.Cmax,
                α = th.α, FB = th.FB, β = th.β, F1 = f1, F2 = f2, t1 = t1, t2 = t2))))
        end
    end
    write_long(joinpath(OUT, "atlas_stein_2019.csv"), rows)
end

# ---------------------------------------------------------------- Lindauer 2017 (mouse model)
let M = load_model("lindauer_2017")
    rows = Tuple{String, Float64, String, Float64}[]
    raw = M.mouse_raw()
    BW = 0.02                                  # kg, Methods (standard mouse 20 g)
    dose_days = [0.0, 5.0, 9.0, 13.0, 17.0]    # test data set design (Methods)
    grid = collect(0:2.0:720.0)                # h, 30 days
    for guard in (false, true)
        sys = M.build_lindauer(; guard)
        for mgkg in (0.0, 0.1, 1.4, 5.0)
            nmol = mgkg * BW * 1e6 / M.MW      # mg/kg * kg -> mg -> nmol (mg * 1e6 / MW)
            prob = ODEProblem(sys, M.problem_pairs(sys, raw), (0.0, last(grid)))
            kw = mgkg > 0 ? (; callback = bolus_callback(sys, sys.A1, dose_days .* 24, nmol)) : (;)
            sol = solve(prob, Rodas5P(); abstol = 1e-14, reltol = 1e-10, saveat = grid, dtmax = 24.0, kw...)
            SciMLBase.successful_retcode(sol) || error("lindauer solve failed $mgkg")
            name = string("mouse_", replace(string(mgkg), "." => "p"), "mgkg", guard ? "_guard" : "")
            CP = sol[sys.A1] ./ raw.V1 .* 1000 .* M.MW ./ 1e6     # nmol / L -> mg/L (s006 $ERROR CP)
            for (i, tt) in enumerate(sol.t)
                push!(rows, (name, tt, "CP", CP[i]))
                push!(rows, (name, tt, "ROb", sol[sys.ROb][i]))
                push!(rows, (name, tt, "ROt", sol[sys.ROt][i]))
                push!(rows, (name, tt, "TV", sol[sys.TV][i]))
            end
        end
    end
    write_long(joinpath(OUT, "atlas_lindauer_2017.csv"), rows)
end
