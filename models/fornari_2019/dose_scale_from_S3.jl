# DIAGNOSTIC ONLY (not used by reproduce()): which dose amount the model LINES of Supplement Figure S3 (labelled 30 mg/kg, dose at t = 24 h)
# imply with the Table S3 PK. Result: 117 umol in Cen (Damage curve: 132), i.e. about 43 mg/kg under the per-kg reading of Fix pk_per_kg,
# so the plotted model line does not correspond to the 20 or 30 mg/kg of the caption (paper issue). The graded conversion is
# 1000/371.25 umol/kg per mg/kg (Fix pk_per_kg, evidence in pk_units_check.jl); DOSE_UMOL_PER_MGKG_S3 is kept for sensitivity.jl.
# The Figure S3 curves are digitised in reference/ancillary/ (blue free platinum in uM, orange Damage). Only one scale factor is
# estimated, by linear least squares of the model curve (linear in dose) on the digitised curve; PK and Damage parameters are as printed.
#   julia --startup-file=no --project=. models/fornari_2019/dose_scale_from_S3.jl
using LinearAlgebra, Printf, Statistics

QC, Vcen, Vper, CL, kdam, krep = 0.06, 0.26, 0.30, 0.42, 0.017, 0.017   # Table S3, Table 1
# states: Cen [umol], Per [umol], Damage ; Damage' = kdam*Cen/Vcen - krep*Damage
A = [-(CL + QC) / Vcen   QC / Vper   0.0;
      QC / Vcen         -QC / Vper   0.0;
      kdam / Vcen        0.0        -krep]
function readcsv(p)
    xs = Float64[]; ys = Float64[]
    for l in eachline(p)
        (startswith(l, "#") || startswith(l, "x,y")) && continue
        a = split(l, ','); push!(xs, parse(Float64, a[1])); push!(ys, parse(Float64, a[2]))
    end
    xs, ys
end
dir = joinpath(@__DIR__, "reference", "ancillary")
td, yd = readcsv(joinpath(dir, "S3_damage.csv"))
tc, yc = readcsv(joinpath(dir, "S3_freePt_uM.csv"))
unit(t) = exp(A * (t - 24.0)) * [1.0, 0.0, 0.0]      # response to 1 umol in Cen at t = 24 h
# damage: least squares for the dose in umol
md = [unit(t)[3] for t in td]
D_umol = dot(md, yd) / dot(md, md)
rd = sqrt(mean((D_umol .* md .- yd) .^ 2)) / 5
# free Pt: log-space least squares on t in [26.5, 47.5] (avoids marker windows at 24-25 h and uses the clean terminal phase), excluding marker columns at 27, 30, 36 h +- 0.4 h
keep = [t >= 26.5 && all(abs(t - m) > 0.4 for m in (27.0, 30.0, 36.0)) for t in tc]
mc = [unit(t)[1] / Vcen for t in tc]
lg = mean(log.(yc[keep]) .- log.(mc[keep]))
C_umol = exp(lg)
@printf "from Damage curve: dose in Cen = %.1f umol for 30 mg/kg  (rms %.4f of axis)\n" D_umol rd
@printf "from free-Pt curve: dose in Cen = %.1f umol for 30 mg/kg\n" C_umol
@printf "ratio Damage/Cp estimates: %.3f\n" D_umol / C_umol
@printf "umol per mg/kg: Damage %.3f, Cp %.3f\n" D_umol / 30 C_umol / 30
@printf "for reference: 0.32 kg and MW 371.25 g/mol gives %.3f umol per mg/kg\n" 0.32 * 1000 / 371.25

println("\n t[h]   Cp_dig   Cp(117)  Cp(132)    Dam_dig  Dam(117) Dam(132)")
for t in (24.1, 24.5, 25.0, 25.97, 27.17, 28.36, 30.76, 35.5, 40.3, 47.5)
    i = argmin(abs.(tc .- t)); j = argmin(abs.(td .- t))
    u = unit(tc[i]); u2 = unit(td[j])
    @printf "%6.2f %8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f\n" tc[i] yc[i] 117u[1]/Vcen 132u[1]/Vcen yd[j] 117u2[3] 132u2[3]
end
