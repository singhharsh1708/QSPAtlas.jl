# Stein et al. 2019, CPT Pharmacometrics Syst Pharmacol 8:285-295 (PMC6539725)
# Cellular kinetic model of tisagenlecleucel (CD19 CAR-T) transgene copies in blood.
#
# Paper form (Eq. 1, "base model"), C = transgene copies per ug genomic DNA, t in days after infusion:
#
#   C(t) = R0 exp(rho t)                                         t <  Tmax
#   C(t) = A exp(-alpha (t - Tmax)) + B exp(-beta (t - Tmax))    t >= Tmax
#   rho = log(foldx) / Tmax,  R0 = Cmax / foldx,  A = Cmax (1 - FB),  B = Cmax FB
#
# With comedication (Eq. 2 / Figure 1b / supplement Monolix code) the exponential phase is
# piecewise: rate rho before T1 (first comedication), F1 rho from T1 to T2, F1 F2 rho from T2
# to Tmax; T1, T2 are clipped to Tmax, and the peak becomes CmaxAdjust (supplement lines
# 125-134) instead of Cmax.
#
# Implemented as the semimechanistic compartment model of Figure 2 / supplement lines 6-31,
# whose analytic solution is Eq. 1 exactly:
#
#   effector cells E:  dE/dt = rho_eff(t) E  (t < Tmax),   dE/dt = -alpha E  (t >= Tmax)
#   memory cells   M:  dM/dt = 0             (t < Tmax),   dM/dt = k E - beta M  (t >= Tmax)
#   C = E + M,  E(0) = R0,  M(0) = 0,  k = FB (alpha - beta)
#
# (supplement line 31: C2 = Cmax k/(alpha - beta) = Cmax FB, so k = FB (alpha - beta)).
# The piecewise switch is `ifelse` on t against the parameter Tmax (and T1, T2). The solver
# steps through the switch; the accuracy against the closed form is checked in reproduce.jl
# (`check_model`, which compares the ODE solution with `closed_form` in protocol.jl).
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Cmax is not a system parameter: it enters only through the initial value E(0) = R0 = Cmax / foldx
@parameters foldx Tmax α FB β F1 F2 T1 T2
@variables E(t) M(t) C(t)

# rho is a derived quantity, not a free parameter (Methods: rho = log(foldx)/Tmax)
ρ = log(foldx) / Tmax
k = FB * (α - β)

@mtkcompile CART = System(
    [D(E) ~ ifelse(t < Tmax, ifelse(t < T1, 1.0, ifelse(t < T2, F1, F1 * F2)) * ρ, -α) * E,
        D(M) ~ ifelse(t < Tmax, 0.0, k * E - β * M),
        C ~ E + M],
    t)

# ---- Table 1 typical values (fixed effects; geometric means of the log-normal parameters) ----
# Table 1 "Fixed effect" rows. Units: foldx [-], Tmax [day], Cmax [copies/ug], F [-], rates [1/day].
const TYPICAL = (
    foldx = 3900.0,    # Table 1: fold expansion, RSE 30 %
    Tmax = 9.3,        # Table 1: day, RSE 4.2 %
    Cmax = 24000.0,    # Table 1: DNA copies/ug, RSE 20 %
    Ftoci = 1.2,       # Table 1: effect of tocilizumab on rho (RSE 7.5 %); NOT used in the base model
    Fster = 1.0,       # Table 1: effect of corticosteroids on rho (RSE 9 %); NOT used in the base model
    α = 0.16,          # Table 1: 1/day (RSE 11 %), half-life ln2/alpha = 4.3 d (text)
    FB = 0.0079,       # Table 1: fraction of Cmax declining with beta (RSE 15 %)
    β = 0.0032,        # Table 1: 1/day (RSE 23 %), half-life ln2/beta = 220 d (text)
)

# ---- Table 1 random effects: omega, read as the STANDARD DEVIATION of the normal eta ---------
# Table 1 "Random effect" rows give one number per parameter, no unit, no "CV %"; all parameters
# are log-normal (Methods "Random-effect model"), Cmax_i = Cmax exp(eta_i) (supplement Eq. 4)
# with eta ~ N(0, omega^2) (supplement line 77: "random effect eta_i with variance omega^2").
# Monolix reports omega as a standard deviation, so the table numbers are omega, not omega^2.
# Correlations between etas are not reported: diagonal Omega (Monolix default).
const OMEGA = (
    foldx = 2.4,       # Table 1, RSE 9.5 %, eta shrinkage 0.39
    Tmax = 0.38,       # Table 1, RSE 7.9 %, shrinkage 0.14
    Cmax = 0.65,       # Table 1, RSE 10 %, shrinkage 0.29
    α = 0.91,          # Table 1, RSE 8.8 %, shrinkage 0.27
    FB = 0.8,          # Table 1, RSE 15 %, shrinkage 0.53
    β = 0.86,          # Table 1, RSE 23 %, shrinkage 0.82
)

# ---- residual error ---------------------------------------------------------------------
# Table 1 "Residual error a = 0.56" (RSE 3.3 %). Methods: "log-transformed data were modeled with
# a proportional error model (constant error model in log-space)": log y = log f + a eps,
# eps ~ N(0, 1). a is read as the standard deviation (Monolix "a" of a constant error model).
const RESIDUAL_SD_LOG = 0.56

# Data below the limit of quantification (BLQ): LLOQ = 10 copies per reaction, about 50
# transgene copies per ug genomic DNA (Methods "Sample analysis"). The fit treated BLQ as
# fixed-point censored data (Methods "Analysis").
const LLOQ = 50.0
