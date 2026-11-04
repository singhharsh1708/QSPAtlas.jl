# Hutchinson et al. 2016, CPT:PSP 5:636-645 (DOI 10.1002/psp4.12142, PMC5192995).
# Logistic tumour growth T(t) with a dynamic vessel-dependent carrying capacity V(t) (Hahnfeldt-type),
# antiangiogenic vessel kill delta_V, a transient "normalisation window" factor N(t), and an optional
# cytotoxic drug C(t) with resistance (Eqs 7-11). Time is in days after tumour cell injection;
# treatment and the first measurement are at day 38 (Methods, Preclinical data).
#
# One system covers all equations of the paper:
#   Eqs 1-2 (control)      delta_V = 0, N = 1                      (Fig 2a,b)
#   Eqs 3-5 (treatment)    N = Nmax in (tnorm1, tnorm2), else 1    (Fig 2c-f, Fig 6a)
#   Eqs 7-11 (+ cytotoxic) delta_T > 0, alphaC_on > 0              (Fig 6b-d)
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
    hutchinson_system(; name)

Compiled MTK system. States `T` and `V` [mm^3]; observed `N` (dimensionless window factor) and `C`
(intratumoural cytotoxic concentration). All parameter defaults are the Figure 2a control values;
callers override them per scenario with `hutchinson_problem`.
"""
function hutchinson_system(; name = :hutchinson)
    @parameters begin
        alpha_T = 0.1        # [1/day] tumour growth rate, Eq. 1; Fig 2 caption (0.1), Table 1 final model 0.109
        alpha_V = 0.09       # [1/day] vessel growth rate, Eq. 2; Fig 2 caption (0.09), Table 1 final model 0.119
        delta_V = 0.0        # [1/day] vessel kill rate, Eq. 4; 0 for the control Eqs 1-2
        Nmax = 1.0           # [-] carrying-capacity enhancement in the window, Eq. 5 / Eq. 10
        tnorm1 = 52.0        # [day] window start, Eq. 5; Fig 2e caption 52, Table 1 53.4
        tnorm2 = 62.0        # [day] window end, Eq. 5; Fig 2e caption 62, Table 1 59.2
        Noff = 0.0           # [-] 0 = Eq. 5 as printed (N = Nmax in window); 1 = supplement code reading N = Nmax + 1 (diagnostic only)
        delta_T = 0.0        # [mL/mg/day] cytotoxic kill rate, Eq. 7; 0.12 in Section "Theoretical administration of chemotherapy"
        lambda = 0.0         # [1/day] resistance rate, Eq. 7; 0.08 in the chemotherapy section
        t_res0 = 0.0         # [day] origin of the resistance clock exp(-lambda (t - t_res0)); Eq. 7 prints exp(-lambda t), i.e. 0
        alphaC_on = 0.0      # [mg/mL/mm^3] alpha~_C, Eq. 11; 1 in the chemotherapy section
        k_C = 0.9            # [1/day] intraday decay of C, Eq. 9; 0.9 in the chemotherapy section
        tCon = 0.0           # [day] chemotherapy on, Eq. 11
        tCoff = 0.0          # [day] chemotherapy off, Eq. 11
    end
    @variables T(t) V(t) N(t) C(t)
    in_window = (t > tnorm1) & (t < tnorm2)
    in_chemo = (t > tCon) & (t < tCoff)
    eqs = [
        N ~ ifelse(in_window, Nmax + Noff, 1.0),                              # Eq. 5 / 10
        C ~ ifelse(in_chemo, alphaC_on, 0.0) * N * V * exp(-k_C * (t - floor(t))),  # Eq. 9 and 11, mod(t,1) = t - floor(t)
        D(T) ~ alpha_T * T * (1 - T / (N * V)) - delta_T * C * T * exp(-lambda * (t - t_res0)),  # Eq. 3 / 7
        D(V) ~ alpha_V * T^(2 / 3) * V^(1 / 3) - delta_V * V,                    # Eq. 2 / 4, beta = 2/3, gamma = 1/3
    ]
    return mtkcompile(System(eqs, t; name))
end

"""
    initial_state(T0, K; K_is_V_over_T = false)

Initial conditions `(T0, V0)`. The paper prints `K = T0/V0` (Methods, before Eq. 6), which gives `V0 = T0/K`.
With `K_is_V_over_T = true`, `V0 = K * T0` (supplement code `V_0 = K*X`, Table 1 "initial ratio of V/T", Fig 2b).
"""
initial_state(T0, K; K_is_V_over_T = false) = (T0, K_is_V_over_T ? K * T0 : T0 / K)
