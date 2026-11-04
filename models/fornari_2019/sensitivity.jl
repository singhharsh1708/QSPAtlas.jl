# Sensitivity / diagnostic runs for fornari_2019 (NOT part of reproduce(); does not change the grade).
# Each row re-scores the frozen reference curves with one reading changed. Nothing here is fitted.
#   julia --startup-file=no --project=. models/fornari_2019/sensitivity.jl
using QSPAtlas, Printf
const m = Module()
Base.include(m, joinpath(@__DIR__, "reproduce.jl"))
refs = load_references(joinpath(@__DIR__, "reference"))
F = m.FIXES
without(ids...) = [f for f in F if !(f.id in ids)]

function score(ref, x, y)
    c = compare_curve(ref, x, y; tolerance = m.TOLERANCE)
    return @sprintf("%.3f%s", c.rms, c.passed ? "*" : " ")
end
function sol_series(s, v)
    sol = m.solve(s.prob, m.Rodas5P(); abstol = 1e-9, reltol = 1e-8, saveat = s.saveat)
    return sol.t ./ 24, sol[getproperty(s.sys, v)]
end

println("Figure 3 (rms in axis-box units, * = pass at 0.015/0.03); columns MPP CMP MEP Neut Plt")
ids3 = ["fig3_MPP", "fig3_CMP", "fig3_MEP", "fig3_Neut", "fig3_Plt"]
vars3 = [:MPP, :CMP, :MEP, :Neut, :Plt]
function row3(label, s)
    cols = [score(refs[i], sol_series(s, v)...) for (i, v) in zip(ids3, vars3)]
    println(rpad(label, 62), join(cols, "  "))
end
row3("main (fixes; 0.862 umol per mg/kg = 0.320 kg, MW 371.25; dose at 24 h)", m.rat_problem(F))
row3("DIAGNOSTIC dose 3.90 umol per mg/kg (Fig S3 PK-implied, ~1.45 kg)", m.rat_problem(F; dose_scale = 3.9))
row3("dose 4.41 umol per mg/kg (Damage-curve estimate of Fig S3)", m.rat_problem(F; dose_scale = 4.41))
row3("first dose at t = 0 instead of 24 h", m.rat_problem(without("dose_time")))
row3("dose 3.00 umol per mg/kg", m.rat_problem(F; dose_scale = 3.0))
row3("dose 3.50 umol per mg/kg", m.rat_problem(F; dose_scale = 3.5))

println("\nFigure 4 (columns: 4a neut, 4a plt, 4b neut, 4b plt)")
ids4 = ["fig4a_neut", "fig4a_plt", "fig4b_neut", "fig4b_plt"]
function human_variant(label; slopes_a = :rat, slopes_b = :table2, kw...)
    cols = String[]
    for (sl, vs, rids) in ((slopes_a, (:Neut, :Plt), ids4[1:2]), (slopes_b, (:Neut, :Plt), ids4[3:4]))
        s = m.human_problem(F; slopes = sl, kw...)
        for (v, rid) in zip(vs, rids)
            push!(cols, score(refs[rid], sol_series(s, v)...))
        end
    end
    println(rpad(label, 62), join(cols, "  "))
end
human_variant("main: 750 mg over 30 min at day 0 and 21, Table 2 slopes")
human_variant("single dose only (S11/S12 wording)"; dose_times_h = [0.0])
human_variant("bolus instead of infusion"; infusion = 0.0)
human_variant("Eq 10 slopes (Slope_R x IC50_R/IC50_H) in (b)"; slopes_b = :eq10)
human_variant("Neut0 4.7 (figure baseline), others midpoint"; phys = m.human_physiology(Neut0 = 4.7))
human_variant("Mono0 0.5"; phys = m.human_physiology(Mono0 = 0.5))
human_variant("Mono0 0.2"; phys = m.human_physiology(Mono0 = 0.2))
human_variant("Ret0 40, RBC0 4100"; phys = m.human_physiology(Ret0 = 40.0, RBC0 = 4100.0))
human_variant("Ret0 115, RBC0 5900"; phys = m.human_physiology(Ret0 = 115.0, RBC0 = 5900.0))
human_variant("exposure x 0.5 (dose 375 mg)"; dose_umol = m.HUMAN_DOSE_UMOL / 2)
human_variant("exposure x 2 (dose 1500 mg)"; dose_umol = m.HUMAN_DOSE_UMOL * 2)
human_variant("AUC-based CL (Calvert CL 150 mL/min = 9 L/h)"; pk = (; QC = 3.46, Vcen = 15.5, Vper = 9.86, CL = 9.0))

# two doses, day 0 and day 21 (the paper's red lines show a second trough near day 35); a reading taken from the figure, so :from_target in kind
human_variant("three doses, day 0, 21, 42"; dose_times_h = [0.0, 504.0, 1008.0])
