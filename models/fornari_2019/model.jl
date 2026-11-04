# Fornari et al. 2019, CPT Pharmacometrics Syst Pharmacol 8:858-868 (PMC6875710).
# Haematopoiesis systems pharmacology model of carboplatin-induced myelosuppression, rat and human.
#
# Time unit: hour. Cell counts: 1e9 cells/L. Free platinum: umol/L (amounts in umol). Damage: umol/L (Slope is 1/(umol/L)).
#
# Equation numbers refer to the main text, "S" numbers to the supplement (PSP4-8-858-s012.pdf), "Table S" to s008-s011,
# "code" to the parameter script in PSP4-8-858-s013.pdf (its ODE right-hand sides are missing from the PDF).
#
# Rate constants are NOT free parameters in the paper: they are derived from baselines, mean transit times, circulating rates and the
# dimensionless lambda ratios by solving the steady state (supplement S6 Eq S3/S4, S7 Eq S7). `derive_rates` does that algebra.

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# ---------------------------------------------------------------------------------------------------------------------
# Physiology-dependent parameters (Table 1) as NamedTuples
# ---------------------------------------------------------------------------------------------------------------------

"""Rat column of Table 1 (and the parameter script in s013)."""
rat_physiology() = (;
    kcirc_Neut = 0.17,      # Table 1, 1/h
    kcirc_Mono = 0.06,      # Table 1
    kcirc_Plt = 0.01,       # Table 1
    kcirc_RBC = 0.0007,     # Table 1
    Neut0 = 1.32,           # Table 1, 1e9 cells/L
    Mono0 = 0.16,           # Table 1
    Plt0 = 967.0,           # Table 1
    Ret0 = 235.0,           # Table 1
    RBC0 = 8161.0,          # Table 1
    MPP0 = 69.0,            # Table 1
    CMP0 = 214.0,           # Table 1
    MEP0 = 156.0,           # Table 1
    MTT_Neut = 61.5,        # Table 1, h
    MTT_Mono = 62.2,        # Table 1
    MTT_Plt = 74.0,         # Table 1
    MTT_Ret = 75.0,         # Table 1
    lambda1 = 2.0,          # Table 1: T2Ret0/T1Ret0
    lambda2 = 2.0,          # Table 1: T2Plt0/T1Plt0
    lambda3 = 1.8,          # Table 1: (ktrNeut+ktrMono)/kprolCMP
    lambda4 = 1.8,          # Table 1: (ktrRet+ktrPlt)/kprolMEP
    lambda5 = 1.8,          # Table 1: (ktrCMP+ktrMEP)/kprolMPP
    gamma_stem = 0.07,      # Table 1 (code: gaMPP)
    gamma_matCMP = 0.6,     # Table 1 (code: gaMPPMatCMP)
    gamma_matMEP = 0.3,     # Table 1 (code: gaMPPMatMEP)
    gamma_prolTrans = 0.7,  # Table 1 (code: gaT1T2)
    alphaWBC = 0.5,         # code (alphaWBC); equals the 1/2 weights printed in Eqs 11-13
    kdam = 0.017,           # Table 1, 1/h
    krep = 0.017,           # Table 1, 1/h
)

"""
Human column of Table 1. Neut0, Mono0, Ret0 and RBC0 are printed only as literature ranges (2-7, 0.2-10, 40-115, 4100-5900); the
keywords take the value to use (default: the midpoint of the printed range, see Fix `human_baselines`).
"""
human_physiology(; Neut0 = 4.5, Mono0 = 5.1, Ret0 = 77.5, RBC0 = 5000.0) = (;
    kcirc_Neut = 0.1,       # Table 1
    kcirc_Mono = 0.04,      # Table 1
    kcirc_Plt = 0.0052,     # Table 1
    kcirc_RBC = 0.00037,    # Table 1
    Neut0, Mono0, Ret0, RBC0,   # Table 1 ranges, midpoints by default
    Plt0 = 345.0,           # Table 1
    MPP0 = 1.3,             # Table 1 (allometric scaling, S9 Eq S8)
    CMP0 = 20.9,            # Table 1 (S9: 16 x MPP0 would give 20.8)
    MEP0 = 15.0,            # Table 1 (S9: 11 x MPP0 would give 14.3)
    MTT_Neut = 210.0,       # Table 1, h
    MTT_Mono = 121.5,       # Table 1
    MTT_Plt = 168.0,        # Table 1
    MTT_Ret = 66.0,         # Table 1
    lambda1 = 2.0, lambda2 = 2.0, lambda3 = 1.8, lambda4 = 1.8, lambda5 = 1.8,   # Table 1 "as in the rat"
    gamma_stem = 0.07, gamma_matCMP = 0.6, gamma_matMEP = 0.3, gamma_prolTrans = 0.7,   # Table 1 "as in the rat"
    alphaWBC = 0.5,         # code
    kdam = 0.017, krep = 0.017,   # Table 1 "as in the rat"
)

# Drug-dependent parameters (Table 2)
rat_drug() = (; Slope_MPP = 2.05, Slope_CMP = 1.47, Slope_MEP = 2.19, delta_Ret = 2.8, delta_Plt = 0.54)   # Table 2 rat, 1/uM; delta dimensionless
human_drug_table2() = (; Slope_MPP = 0.79, Slope_CMP = 0.57, Slope_MEP = 0.66, delta_Ret = 2.8, delta_Plt = 0.54)   # Table 2 human column "IC50 scaled value"
# Eq 10 applied to the rat slopes with Table S2 IC50 (CMP: 2.83 rat, 4.04 human; MEP: 0.86 rat, 1.56 human). Table S2 has no MPP IC50;
# the MPP slope is scaled with the CMP ratio here (assumption, used only in the sensitivity analysis, never in the graded run).
human_drug_eq10() = (; Slope_MPP = 2.05 * 2.83 / 4.04, Slope_CMP = 1.47 * 2.83 / 4.04, Slope_MEP = 2.19 * 0.86 / 1.56, delta_Ret = 2.8, delta_Plt = 0.54)

# PK parameters
rat_pk() = (; QC = 0.06, Vcen = 0.26, Vper = 0.30, CL = 0.42)   # Table S3 (code lists L/h, L, L, L/h); used per kg (L/h/kg, L/kg) with the dose in umol/kg, Fix pk_per_kg
# Human: NOT in the paper. S11 cites Zandvliet et al. Br J Clin Pharmacol 2008;66:485 (PMC2561111, Table 2): CL = 0.76 CrCl + 1.5 L/h,
# V1 15.5 L, Q 3.46 L/h, V2 9.86 L. S11 gives GFR = 125 mL/min = 7.5 L/h, assumed equal to CrCl in that formula.
human_pk() = (; QC = 3.46, Vcen = 15.5, Vper = 9.86, CL = 0.76 * 7.5 + 1.5)

# ---------------------------------------------------------------------------------------------------------------------
# Steady-state algebra (S6 Eq S3, S4; S7 Eq S7)
# ---------------------------------------------------------------------------------------------------------------------
function derive_rates(p)
    a_Neut = 3 / p.MTT_Neut;  a_Mono = 3 / p.MTT_Mono       # Eq 6
    a_Ret = 3 / p.MTT_Ret;    a_Plt = 3 / p.MTT_Plt         # Eq 8
    kcirc_Ret = p.kcirc_RBC * p.RBC0 / p.Ret0               # S4 B
    # transit baselines (S3): T2 = T3 = kcirc*cells/a; T1 = T2/lambda for the proliferative Ret and Plt chains
    T_Neut = p.kcirc_Neut * p.Neut0 / a_Neut
    T_Mono = p.kcirc_Mono * p.Mono0 / a_Mono
    T_Ret = kcirc_Ret * p.Ret0 / a_Ret
    T_Plt = p.kcirc_Plt * p.Plt0 / a_Plt
    T1_Ret = T_Ret / p.lambda1
    T1_Plt = T_Plt / p.lambda2
    ktr_Neut = p.kcirc_Neut * p.Neut0 / p.CMP0              # S4
    ktr_Mono = p.kcirc_Mono * p.Mono0 / p.CMP0
    kprol_Ret = a_Ret * (1 - 1 / p.lambda1)                 # S4
    kprol_Plt = a_Plt * (1 - 1 / p.lambda2)
    # T1 steady state: kprol T1 + ktr MEP0 - a T1 = 0  (S4: ktr = a T1/MEP0 * T1/T)
    ktr_Ret = a_Ret * T1_Ret / (p.lambda1 * p.MEP0)
    ktr_Plt = a_Plt * T1_Plt / (p.lambda2 * p.MEP0)
    kprol_CMP = (ktr_Neut + ktr_Mono) / p.lambda3           # S7 Eq S7 (lambda3)
    kprol_MEP = (ktr_Ret + ktr_Plt) / p.lambda4             # lambda4
    ktr_CMP = (ktr_Neut + ktr_Mono - kprol_CMP) * p.CMP0 / p.MPP0   # S4
    ktr_MEP = (ktr_Ret + ktr_Plt - kprol_MEP) * p.MEP0 / p.MPP0
    kprol_MPP = (ktr_CMP + ktr_MEP) / p.lambda5             # lambda5
    kstem = (ktr_CMP + ktr_MEP - kprol_MPP) * p.MPP0        # S4
    return (; a_Neut, a_Mono, a_Ret, a_Plt, kcirc_Ret, T_Neut, T_Mono, T_Ret, T_Plt, T1_Ret, T1_Plt, ktr_Neut, ktr_Mono, ktr_Ret, ktr_Plt,
        kprol_Ret, kprol_Plt, kprol_CMP, kprol_MEP, kprol_MPP, ktr_CMP, ktr_MEP, kstem)
end

# ---------------------------------------------------------------------------------------------------------------------
# The system
# ---------------------------------------------------------------------------------------------------------------------
"""
    haematopoiesis_system(phys, drug, pk; structure = :corrected, name = :haem)

MTK `System` for Eqs 1-9, 11-13 with the carboplatin PK of S5 (`Cen`, `Per` in umol; dose enters `Cen` by callback, or through the infusion
rate parameter `rin` [umol/h]). Returns `(sys, u0)` where `u0` is the steady-state initial condition (Damage = 0, PK empty).

`structure = :as_printed` throws: the printed equations are not closed (Eq 1 uses an undefined `fdbk_mat_MPP` and has no MPP factor on
the maturation terms, so the steady state of S6 does not hold; Eq 9 has no RBC equation; Eq 7 platelet T1 uses a_Ret). `:corrected`
applies the Fixes `eq1_maturation`, `rbc_equation` and `plt_transit_rate` described in protocol.jl / report.md.
"""
function haematopoiesis_system(phys, drug, pk; structure::Symbol = :corrected, name = :haem)
    structure === :corrected || error("model as printed cannot be run: Eq 1 contains fdbk_mat_MPP, which is defined nowhere (Eq 12 defines fdbk_mat_CMP and fdbk_mat_MEP), " *
        "and its maturation terms lack the factor MPP, so the steady state of supplement S6 is violated; Eq 9 omits the RBC equation; Eq 7 uses a_Ret in the platelet T1 equation. " *
        "The supplementary ODE code (s013) is truncated before the right-hand sides.")
    r = derive_rates(phys)
    @variables MPP(t) CMP(t) MEP(t) N1(t) N2(t) N3(t) M1(t) M2(t) M3(t) R1(t) R2(t) R3(t) P1(t) P2(t) P3(t)
    @variables Neut(t) Mono(t) Ret(t) Plt(t) RBC(t) Damage(t) Cen(t) Per(t) Cp(t)
    @parameters rin = 0.0   # infusion rate into Cen, umol/h (used for the human 30 min infusion; zero otherwise)
    α = phys.alphaWBC
    # feedback functions, Eqs 11-13
    fdbk_stem = (α * phys.CMP0 / CMP + (1 - α) * phys.MEP0 / MEP)^phys.gamma_stem
    fdbk_mat_CMP = (α * phys.Neut0 / Neut + (1 - α) * phys.Mono0 / Mono)^phys.gamma_matCMP
    fdbk_mat_MEP = (phys.RBC0 / RBC)^phys.gamma_matMEP
    fdbk_prol_Ret = (α * phys.Ret0 / Ret + (1 - α) * phys.RBC0 / RBC)^phys.gamma_prolTrans
    fdbk_prol_Plt = (phys.Plt0 / Plt)^phys.gamma_prolTrans
    ke = pk.CL / pk.Vcen; k12 = pk.QC / pk.Vcen; k21 = pk.QC / pk.Vper
    eqs = [
        Cp ~ Cen / pk.Vcen,                                                          # S5
        D(Cen) ~ -ke * Cen - k12 * Cen + k21 * Per + rin,                            # Eq S1
        D(Per) ~ k12 * Cen - k21 * Per,
        D(Damage) ~ phys.kdam * Cp - phys.krep * Damage,                             # Eq 3
        # Eq 1 (corrected: maturation terms proportional to MPP with the Eq 4 feedbacks)
        D(MPP) ~ r.kstem * fdbk_stem + r.kprol_MPP * (1 - drug.Slope_MPP * Damage) * MPP -
                 r.ktr_CMP * fdbk_mat_CMP * MPP - r.ktr_MEP * fdbk_mat_MEP * MPP,
        # Eq 4
        D(CMP) ~ r.kprol_CMP * (1 - drug.Slope_CMP * Damage) * CMP + r.ktr_CMP * fdbk_mat_CMP * MPP - r.ktr_Neut * CMP - r.ktr_Mono * CMP,
        D(MEP) ~ r.kprol_MEP * (1 - drug.Slope_MEP * Damage) * MEP + r.ktr_MEP * fdbk_mat_MEP * MPP - r.ktr_Ret * MEP - r.ktr_Plt * MEP,
        # Eq 5, 6
        D(N1) ~ r.ktr_Neut * CMP - r.a_Neut * N1, D(N2) ~ r.a_Neut * N1 - r.a_Neut * N2, D(N3) ~ r.a_Neut * N2 - r.a_Neut * N3,
        D(M1) ~ r.ktr_Mono * CMP - r.a_Mono * M1, D(M2) ~ r.a_Mono * M1 - r.a_Mono * M2, D(M3) ~ r.a_Mono * M2 - r.a_Mono * M3,
        # Eq 7, 8 (proliferative T1, T2; non-proliferative T3)
        D(R1) ~ r.kprol_Ret * fdbk_prol_Ret * R1 - drug.delta_Ret * drug.Slope_MEP * r.kprol_Ret * Damage * R1 + r.ktr_Ret * MEP - r.a_Ret * R1,
        D(R2) ~ r.kprol_Ret * fdbk_prol_Ret * R2 - drug.delta_Ret * drug.Slope_MEP * r.kprol_Ret * Damage * R2 + r.a_Ret * R1 - r.a_Ret * R2,
        D(R3) ~ r.a_Ret * R2 - r.a_Ret * R3,
        D(P1) ~ r.kprol_Plt * fdbk_prol_Plt * P1 - drug.delta_Plt * drug.Slope_MEP * r.kprol_Plt * Damage * P1 + r.ktr_Plt * MEP - r.a_Plt * P1,
        D(P2) ~ r.kprol_Plt * fdbk_prol_Plt * P2 - drug.delta_Plt * drug.Slope_MEP * r.kprol_Plt * Damage * P2 + r.a_Plt * P1 - r.a_Plt * P2,
        D(P3) ~ r.a_Plt * P2 - r.a_Plt * P3,
        # Eq 9 (+ RBC equation, S4: kcirc_Ret = kcirc_RBC RBC0/Ret0)
        D(Neut) ~ r.a_Neut * N3 - phys.kcirc_Neut * Neut,
        D(Mono) ~ r.a_Mono * M3 - phys.kcirc_Mono * Mono,
        D(Ret) ~ r.a_Ret * R3 - r.kcirc_Ret * Ret,
        D(Plt) ~ r.a_Plt * P3 - phys.kcirc_Plt * Plt,
        D(RBC) ~ r.kcirc_Ret * Ret - phys.kcirc_RBC * RBC,
    ]
    sys = mtkcompile(System(eqs, t; name))
    u0 = [sys.MPP => phys.MPP0, sys.CMP => phys.CMP0, sys.MEP => phys.MEP0,
        sys.N1 => r.T_Neut, sys.N2 => r.T_Neut, sys.N3 => r.T_Neut, sys.M1 => r.T_Mono, sys.M2 => r.T_Mono, sys.M3 => r.T_Mono,
        sys.R1 => r.T1_Ret, sys.R2 => r.T_Ret, sys.R3 => r.T_Ret, sys.P1 => r.T1_Plt, sys.P2 => r.T_Plt, sys.P3 => r.T_Plt,
        sys.Neut => phys.Neut0, sys.Mono => phys.Mono0, sys.Ret => phys.Ret0, sys.Plt => phys.Plt0, sys.RBC => phys.RBC0,
        sys.Damage => 0.0, sys.Cen => 0.0, sys.Per => 0.0]
    return sys, u0
end
