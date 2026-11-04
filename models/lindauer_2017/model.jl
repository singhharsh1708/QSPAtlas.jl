# Lindauer et al. 2017 (CPT Pharmacometrics Syst Pharmacol 6:11-20, PMC5270293):
# translational PK / PD-1 receptor occupancy / tumour growth inhibition model of pembrolizumab
# (mouse surrogate DX400 -> human). Physiologically based tumour tissue compartment of Shah & Betts
# coupled to a two-compartment plasma PK model with target-mediated binding, PD-1 upregulation and a
# Simeoni-type tumour growth model with the drug effect driven by receptor occupancy in the tumour.
#
# Equations: supplement s001 ("Mathematical representation of the model") and the NONMEM control
# stream s006 (mouse model, ADVAN13). Where the two differ, the control stream is followed and the
# difference is listed in report.md. Units inside the ODE system follow the control stream:
# time [h], amounts [nmol], concentrations [nM], volumes [L], tumour volume TV [mm^3 = uL].
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

const AV = 6.0221415e23   # s006 line "AV = 6.0221415E23" (Avogadro)
const MW = 149000.0       # s006 "MW = 149000 g/mol" (antibody)
const PSI = 20.0          # s006 "PSI = 20" (Simeoni), supplement tumour volume equation
const F_ES = 0.005        # V_es = 0.5 % of V_tot, Table 1
const F_IS = 0.55         # V_is = 55 % of V_tot, Table 1
const F_VS = 0.07         # V_vs = 7 % of V_tot, Table 1

"""
Raw (paper-unit) parameter set. Field units: V1, V2 [mL]; Q, CL [mL/day]; Vmax [ug/day]; KM [ug/mL];
PLQ [L/h/L]; Lfrac [-] (lymph flow as fraction of plasma flow); CLup [L/h/L]; Kdeg [1/h]; FcRni [uM];
KonFcRn, KonPD1 [1e6/M/h]; KoffFcRn, KoffPD1, KdegPD1 [1/h]; N_Tcell [1/uL]; V_blood [uL];
EC50TP [nM]; L0 [1/day]; L1 [mm^3/day] (Inf = exponential growth only); W0 [mm^3];
SLtg [1/(% day)] (drug effect DE = SLtg * RO^gamma per day, RO in %).
"""
const RAW_FIELDS = (:V1, :V2, :Q, :Vmax, :KM, :CL, :PLQ, :Lfrac, :CLup, :Kdeg, :vref, :vref_is, :FcRni, :FR,
    :KonFcRn, :KoffFcRn, :N_Tcell, :V_blood, :N_PD1_TC, :Tmulti, :EMAXTP, :EC50TP, :KonPD1, :KoffPD1, :K_IVIV,
    :KdegPD1, :L0, :L1, :W0, :SLtg, :gamma)

# Parameters shared by mouse and man ("constant across species" in Table 1). Estimated ones from Table S2.
const SPECIES_INVARIANT = (
    PLQ = 12.7,        # Table 1, Shah et al.
    Lfrac = 0.002,     # Table 1: lymph flow 0.2 % of plasma flow
    CLup = 0.0366,     # Table 1
    Kdeg = 42.9,       # Table 1
    vref = 0.842,      # Table 1
    vref_is = 0.2,     # Table 1
    FcRni = 49.8,      # Table 1 [uM]
    FR = 0.715,        # Table 1
    N_Tcell = 1000.0,  # Table 1 [1/uL]
    N_PD1_TC = 10000.0, # Table 1 (assumed)
    Tmulti = 4.32,     # Table S2 estimate (Table 1: 4.3)
    EMAXTP = 94.7,     # Table S2 / Table 1
    EC50TP = 1.46,     # Table S2 / Table 1 [nM]
    K_IVIV = 1.0,      # s006 THETA(25)
    gamma = 2.28,      # Table S2 estimate (Table 1: 2.3)
)

"Mouse parameter set, Table 1 'Value in mouse' / Table S2 estimates."
mouse_raw() = (; SPECIES_INVARIANT...,
    V1 = 1.26, V2 = 0.819, Q = 4.82, Vmax = 0.518, KM = 0.366, CL = 0.334,   # Table 1, Table S2 (mouse PK)
    KonFcRn = 80.6, KoffFcRn = 6.55,                                         # Table 1
    V_blood = 1400.0,                                                        # s006 THETA(18) [uL]
    KonPD1 = 340.0, KoffPD1 = 0.106, KdegPD1 = 0.0194,                       # Table 1 (s006: 340.2, 0.10584)
    L0 = 0.113, L1 = 187.0, W0 = 170.0, SLtg = 1.98e-5,                      # Table 1 / Table S2
)

# Table S3: human melanoma scenarios (baseline volume 41.5 mL = 41500 mm^3; L1 not available, exponential growth)
const SCENARIOS = Dict(
    "fast_allom"   => (L0 = 0.0088, SLtg = 2.575e-6),   # Table S3 row 1
    "fast_growth"  => (L0 = 0.0088, SLtg = 1.542e-6),   # row 2
    "med_allom"    => (L0 = 0.0036, SLtg = 2.575e-6),   # row 3
    "med_growth"   => (L0 = 0.0036, SLtg = 6.247e-7),   # row 4
    "slow_allom"   => (L0 = 0.0017, SLtg = 2.575e-6),   # row 5
    "slow_growth"  => (L0 = 0.0017, SLtg = 2.886e-7),   # row 6
)

"""
Human parameter set. `scenario` is a key of `SCENARIOS`; `V_blood` [uL] is not given for man anywhere in the
paper (see report.md); the default is the only printed value (mouse, s006).
"""
human_raw(; scenario = "fast_allom", V_blood = 1400.0, KdegPD1 = 0.00246) = (; SPECIES_INVARIANT...,
    V1 = 2877.0, V2 = 2854.0, Q = 384.0, Vmax = 114.0, KM = 0.078, CL = 167.0,   # Table 1 'Value in man' (actual human PK)
    KonFcRn = 792.0, KoffFcRn = 23.9,                                            # Table 1
    V_blood,                                                                     # not printed for man
    KonPD1 = 2880.0, KoffPD1 = 0.144, KdegPD1,                                   # Table 1 (in vitro pembrolizumab; allometric Kdeg)
    L0 = SCENARIOS[scenario].L0, L1 = Inf, W0 = 41500.0, SLtg = SCENARIOS[scenario].SLtg,  # Table S3
)

"Convert a raw parameter set to the ODE-level parameter values and initial conditions (control stream PK block)."
function derive(r)
    V1 = r.V1 / 1000                                  # s006: V1 = TVV1/1000 [L]
    V2 = r.V2 / 1000
    Q = r.Q / (1000 * 24)                             # mL/day -> L/h
    VMAX = r.Vmax / (MW * 24) * 1e3                   # ug/day -> nmol/h
    KM = r.KM / MW * 1e6                              # ug/mL -> nM
    CL = r.CL / (1000 * 24)
    PLQraw = r.PLQ
    Lraw = PLQraw * r.Lfrac                           # Lraw = PLQraw*TVL
    M_PD1_b = r.N_PD1_TC * r.N_Tcell * r.V_blood / AV * 1e9   # nmol
    C_PD1_b = M_PD1_b / V1                            # nM
    kon_pd1 = r.KonPD1 * 1e6 / 1e9 / r.K_IVIV         # 1/nM/h, in vivo kon = in vitro kon / K_IVIV
    V_is_in = F_IS * r.W0 / 1e6
    M_PD1_ti = C_PD1_b * r.Tmulti * V_is_in           # nmol
    kout = r.KdegPD1
    p = (V1 = V1, V2 = V2, Q = Q, VMAX = VMAX, KM = KM, CL = CL, PLQraw = PLQraw, Lraw = Lraw,
        CLupraw = r.CLup, Kdeg = r.Kdeg, vref = r.vref, vrefis = r.vref_is, FR = r.FR,
        KonF = r.KonFcRn * 1e6 / 1e9, KoffF = r.KoffFcRn, CPD1b = C_PD1_b,
        KonP = kon_pd1, KoffP = r.KoffPD1, kdp = r.KdegPD1, kin = M_PD1_ti * kout, kout = kout,
        EMAX = r.EMAXTP, EC50 = r.EC50TP, L0 = r.L0 / 24, L0overL1 = isinf(r.L1) ? 0.0 : r.L0 / r.L1,   # L0/L1 [1/mm^3] formed in day units (unit-free ratio of time units)
       
        SL = r.SLtg / 24, gam = r.gamma)
    u0 = (A1 = 0.0, A2 = 0.0, Avs = 0.0, Ceub = 0.0, Ceb = 0.0, FcRn = r.FcRni * 1000, Ais = 0.0,
        PD1t = 0.0, PD1b = 0.0, Mt = M_PD1_ti, TV = r.W0)
    return (; p, u0)
end

@parameters V1 V2 Q VMAX KM CL PLQraw Lraw CLupraw Kdeg vref vrefis FR KonF KoffF CPD1b KonP KoffP kdp kin kout EMAX EC50 L0 L0overL1 SL gam
@variables A1(t) A2(t) Avs(t) Ceub(t) Ceb(t) FcRn(t) Ais(t) PD1t(t) PD1b(t) Mt(t) TV(t)
@variables Vtot(t) Ves(t) Vis(t) Vvs(t) PLQ(t) Lf(t) CLup(t) C1(t) C2(t) Cvs(t) Cis(t) CPD1t(t) ROt(t) ROb(t) DE(t)

"""
Build the compiled ModelingToolkit system (same structure for mouse and man). Equations follow s006 DES block:
A1/A2 plasma amounts, Avs tumour vascular amount, Ceub/Ceb/FcRn endosomal concentrations, Ais interstitial
amount, PD1t/PD1b drug-receptor complex concentrations (tumour interstitium / blood), Mt tumour PD-1 amount, TV tumour volume.
"""
function build_lindauer(; guard = false)
    eqs = [
        Vtot ~ TV / 1e6, Ves ~ F_ES * Vtot, Vis ~ F_IS * Vtot, Vvs ~ F_VS * Vtot,    # $DES: V_tot, V_es, V_is, V_vs
        PLQ ~ PLQraw * Vvs, Lf ~ Lraw * Vvs, CLup ~ CLupraw * Ves,
        C1 ~ A1 / V1, C2 ~ A2 / V2, Cvs ~ Avs / Vvs, Cis ~ Ais / Vis, CPD1t ~ Mt / Vis,
        ROt ~ 100 * PD1t / CPD1t, ROb ~ 100 * PD1b / CPD1b,
        # supplement: DE = SLtg*RO_tumor^gamma. The control stream (s006) uses SLtg*RO for RO <= 0.01 %; `guard = true` switches that on.
        DE ~ guard ? SL * ifelse(ROt > 0.01, max(ROt, 0.0)^gam, ROt) : SL * max(ROt, 0.0)^gam,
        # central compartment (s006 DADT(1))
        D(A1) ~ -CL * C1 - VMAX * C1 / (KM + C1) - PLQ * C1 + PLQ * Cvs - Q * C1 + Q * C2 -
                KonP * A1 * (CPD1b - PD1b) + KoffP * PD1b * V1,
        D(A2) ~ Q * C1 - Q * C2,                                                       # DADT(2)
        D(Avs) ~ PLQ * C1 - (PLQ - Lf) * Cvs - (1 - vref) * Lf * Cvs - CLup * Cvs + CLup * FR * Ceb,   # DADT(3)
        D(Ceub) ~ CLup / Ves * (Cvs + Cis) - KonF * Ceub * FcRn + KoffF * Ceb - Kdeg * Ceub,           # DADT(4)
        D(Ceb) ~ KonF * Ceub * FcRn - KoffF * Ceb - CLup / Ves * Ceb,                  # DADT(5)
        D(FcRn) ~ -KonF * Ceub * FcRn + KoffF * Ceb + CLup / Ves * Ceb,                # DADT(6)
        D(Ais) ~ (1 - vref) * Lf * Cvs - (1 - vrefis) * Lf * Cis - CLup * Cis + CLup * (1 - FR) * Ceb -
                 KonP * Cis * Vis * (CPD1t - PD1t) + KoffP * PD1t * Vis,               # DADT(7)
        D(PD1t) ~ KonP * Cis * (CPD1t - PD1t) - KoffP * PD1t - kdp * PD1t,             # DADT(8)
        D(PD1b) ~ KonP * C1 * (CPD1b - PD1b) - KoffP * PD1b - kdp * PD1b,              # DADT(9)
        D(Mt) ~ kin * (1 + EMAX * PD1t / (EC50 + PD1t)) - kout * Mt,                   # DADT(10)
        D(TV) ~ L0 * TV / (1 + (L0overL1 * TV)^PSI)^(1 / PSI) - DE * TV,             # DADT(11): s006 L0*A/(1+(L0/L1*A)**PSI)**(1/PSI)
    ]
    @named lindauer = System(eqs, t)
    return mtkcompile(lindauer)
end

"Problem pair lists for `ODEProblem(sys, [pairs], tspan)` from a raw parameter set."
function problem_pairs(sys, r)
    d = derive(r)
    ps = [getproperty(sys, k) => v for (k, v) in pairs(d.p)]
    us = [getproperty(sys, k) => v for (k, v) in pairs(d.u0)]
    return vcat(us, ps)
end
