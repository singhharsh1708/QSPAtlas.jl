# Evidence for Fix `pk_per_kg` (Table S3 volumes and clearances read per kg). Not part of reproduce(); uses no target figure.
#   julia --startup-file=no --project=. models/fornari_2019/pk_units_check.jl
#
# Two readings of Table S3 (QC 0.06, Vcen 0.26, Vper 0.30, CL 0.42):
#   per kg  : L/kg and L/h/kg, dose in umol/kg  -> 30 mg/kg = 30/371.25*1000 = 80.8 umol/kg into Vcen 0.26 L/kg
#   absolute: L and L/h for a 0.320 kg rat     -> 30 mg/kg = 25.9 umol into Vcen 0.26 L
# (1) physiology / cited source, (2) the observed in-house 30 mg/kg free-Pt points of Supplement Fig S3 (non-target figure),
# read by eye from the PDF rendered at 100 dpi (log axis: 1e3 at row 45, 1e-2 at row 515, 94 px per decade; +-3 px = +-0.03 decade).
using LinearAlgebra, Printf, Statistics

QC, Vcen, Vper, CL = 0.06, 0.26, 0.30, 0.42
MW = 371.25
A = [-(CL + QC) / Vcen  QC / Vper; QC / Vcen  -QC / Vper]
cp(t, umol) = (exp(A * t) * [umol, 0.0])[1] / Vcen          # t in h after the bolus

println("(1) Physiology, 0.320 kg rat; Siddik 1987 (the paper's PK data source, Biochem Pharmacol 36:1925, abstract):")
println("    carboplatin plasma clearance 10.3 mL/min/kg, inulin 10.1, renal 9.3: glomerular filtration only.")
for (lab, bw) in (("per kg  ", 1.0), ("absolute", 0.320))
    @printf "    %s: CL %.1f mL/min/kg, Vcen %.2f L/kg, Vss %.2f L/kg\n" lab CL / bw * 1000 / 60 Vcen / bw (Vcen + Vper) / bw
end
println("    (rat total body water about 0.67 L/kg, extracellular water about 0.27 L/kg)")

# in-house 30 mg/kg 'o' points: (time after dose [h], rows in px)
obs = [(1.0, [160, 165, 170]), (3.0, [252, 258, 268]), (6.0, [305, 325, 327]), (12.0, [323, 330, 365]), (24.0, [425, 437])]
val(row) = 10^(3 - (row - 45) / 94)
println("\n(2) Fig S3 in-house 30 mg/kg points vs model (uM)")
println("  t[h]  observed            per kg   absolute")
res = Dict("per kg" => Float64[], "absolute" => Float64[])
for (tt, rows) in obs
    ys = val.(rows)
    pk, ab = cp(tt, 30 / MW * 1000), cp(tt, 30 * 0.320 / MW * 1000)
    @printf "%5.1f  %-18s %7.3f  %7.3f\n" tt join([@sprintf("%.2f", y) for y in ys], " ") pk ab
    append!(res["per kg"], log10.(ys) .- log10(pk)); append!(res["absolute"], log10.(ys) .- log10(ab))
end
for k in ("per kg", "absolute")
    @printf "  %-8s mean log10(obs/model) %+.2f, rms %.2f decade (factor %.2f)\n" k mean(res[k]) sqrt(mean(res[k] .^ 2)) 10^mean(res[k])
end
@printf "\nFig S3 model line implies 117 umol in Cen: per kg that is %.1f mg/kg (caption: 20 and 30 mg/kg)\n" 117 * MW / 1000
