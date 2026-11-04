# Ouerdani et al. 2015, CPT Pharmacometrics Syst Pharmacol 4:660-668 (PMC4716582)
# Tumour growth / angiogenesis inhibition model for pazopanib, as printed.
#
#   dP/dt = lambda * P * (1 - P/K) - alpha * exp(-delta * t) * P      (Eq. 1, tumour size P)
#   dK/dt = b * P^n - gamma * K                                       (Eq. 1, carrying capacity K)
#   alpha = alpha0 * AUC^beta_alpha                                   (Eq. 2)
#   gamma = gamma0 * AUC^beta_gamma                                   (Eq. 3)
#
# P: tumour volume [mm^3] (mice) or sum of longest diameters SLD [mm] (patients).
# K: tumour carrying capacity, same unit as P. t in days (NONMEM $DES uses T in days,
# Table 1 gives all rate constants in 1/day). Treatment is continuous: the drug enters only
# through the constant exposure AUC and the decaying cytotoxic term exp(-delta*t); there is
# no dosing event. t = 0 is the start of treatment.
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters λ b γ0 βγ α0 βα δ n AUC
@variables P(t) K(t)

@mtkcompile TUMOUR = System(
    [D(P) ~ λ * P * (1 - P / K) - (α0 * AUC^βα) * exp(-δ * t) * P,
     D(K) ~ b * P^n - (γ0 * AUC^βγ) * K],
    t)

# ---- Table 1 population estimates (typical values) --------------------------------------
# Preclinical (CAKI-2 mice). n = 1: text "best results were obtained with an n value of 1";
# NONMEM s004.txt $DES: DADT(2) = B*A(1) - GAMMA*A(2).
const PRECLINICAL = (
    λ = 0.166,    # 1/day, Table 1
    K0 = 543.0,   # mm^3,  Table 1
    b = 0.0183,   # 1/day, Table 1
    γ0 = 0.007,   # 1/day, Table 1 (gamma at AUC = 1 ug*h/mL)
    βγ = 0.332,   # -,     Table 1
    α0 = 0.251,   # 1/day, Table 1
    βα = 0.0,     # -,     Table 1 "Fixed to 0" (text: initially 0.0002, then fixed to 0)
    δ = 0.196,    # 1/day, Table 1
    n = 1.0,      # -,     Results text ("best results ... n value of 1"), s004.txt line DADT(2)
)
# AUC [ug*h/mL] of the three dose groups, Eq. 2-3 text: 220.2 (10 mg/kg), 656.8 (30), 1140.8 (100).
const MOUSE_AUC = (dose10 = 220.2, dose30 = 656.8, dose100 = 1140.8)

# Clinical (RCC patients, SLD in mm). n = 0.5: text "n ... was set in this case to the value 0.5";
# NONMEM s001.txt $DES: DADT(2) = B*A(1)**0.5 - GAMMA*A(2).
const CLINICAL = (
    λ = 0.0021,   # 1/day, Table 1
    K0 = 329.0,   # mm,    Table 1
    b = 0.0392,   # Table 1 says 1/day; with n = 0.5 and P, K in mm it must be mm^0.5/day (paper's own unit slip)
    γ0 = 0.0023,  # 1/day, Table 1
    βγ = 0.142,   # -,     Table 1
    α0 = 0.0032,  # 1/day, Table 1
    βα = 0.125,   # -,     Table 1
    δ = 0.0153,   # 1/day, Table 1
    n = 0.5,      # -,     Results text, s001.txt $DES
)
# Mean exposure of the 47 patients from the Emax model: 771.6 ug*h/mL (Results text).
const PATIENT_AUC = 771.6

"""
    tumour_problem(par, P0_value, AUC_value, tspan)

ODEProblem for population parameter set `par` (`PRECLINICAL` or `CLINICAL`), initial tumour
size `P0_value` (not given in the paper, see protocol.jl) and exposure `AUC_value`; K starts at `par.K0`.
"""
function tumour_problem(par, P0_value, AUC_value, tspan)
    return ODEProblem(TUMOUR,
        [P => P0_value, K => par.K0,
         λ => par.λ, b => par.b, γ0 => par.γ0, βγ => par.βγ, α0 => par.α0, βα => par.βα,
         δ => par.δ, n => par.n, AUC => AUC_value],
        tspan)
end
