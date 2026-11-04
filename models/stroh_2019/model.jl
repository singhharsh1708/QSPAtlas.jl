# Stroh et al. 2019, PROBODY (Pb-Tx) QSP model, rebuilt from the paper (Fig 2, Tables 1-2, Methods)
# and checked against the KroneckerBio matrix code of the supplement (PSP4-8-676-s003.pdf).
#
# State vector (concentrations, nM; time in days inside the ODE, rate constants kept in 1/s as printed):
#   3 compartments (central = plasma, peripheral, tumour). Each holds
#     R                 free target
#     6 unbound forms   arm states c (mask closed), o (mask open = "breathing"), m (mask cleaved off):
#                       c_c, c_o, o_o, c_m, o_m, m_m   (c_c = intact, c_m/o_m = mono-cleaved, m_m = parental mAb)
#     9 target-bound    c_oR, o_oR, oR_oR, c_mR, o_mR, oR_m, oR_mR, m_mR, mR_mR   (R after an arm = bound)
#     pb_tmdd_m, pb_tmdd_o (periph/tumour: internalised pools, cleaved / intact) and pb_el (central, periph)
#   Reaction rules (read from the supplement matrices A1, A2 and verified entry by entry in
#   supplement_check/check_matrices.jl):
#     per arm: c -> o  at kopen, o -> c at kclose, c/o/oR -> m/m/mR at kcleave (compartment specific)
#     per arm: o or m + R -> oR or mR at kon1 (first arm) or kon2 (second arm); oR/mR -> o/m + R at koff1
#     complexes are internalised at kendo (to tmdd_m if a bound arm is cleaved, else tmdd_o; supplement A1 rows 34/35, 53/54), free R decays at kendo
#     unbound forms: eliminated at kel (central, peripheral) and exchanged central<->peripheral, central<->tumour
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

const COMPS = (:central, :peripheral, :tumor)
const UNBOUND = (:c_c, :c_o, :o_o, :c_m, :o_m, :m_m)
const BOUND = (:c_oR, :o_oR, :oR_oR, :c_mR, :o_mR, :oR_m, :oR_mR, :m_mR, :mR_mR)
const SPEC = (UNBOUND..., BOUND...)

arms(s::Symbol) = Tuple(Symbol.(split(String(s), "_")))
function canon(a, b)   # unordered arm pair -> species name
    for s in SPEC
        arms(s) == (a, b) && return s
        arms(s) == (b, a) && return s
    end
    error("no species with arms ($a, $b)")
end
isbound(s) = any(a -> a in (:oR, :mR), arms(s))
hasmR(s) = any(a -> a == :mR, arms(s))   # supplement: internalised pool is tmdd_m iff a bound arm is cleaved (oR_m goes to tmdd_o)

"""
    Rxn

One reaction in one compartment or between compartments: `from` (1 or 2 state names incl. :R),
`to` (state names), `rate` (key into the parameter dictionary), `mult` (stoichiometric multiplicity of the
rate), `vscale` (volume ratio `(V_from, V_to)` for exchange reactions, `nothing` otherwise).
"""
struct Rxn
    from::Vector{Symbol}
    to::Vector{Symbol}
    rate::Symbol
    mult::Int
    vscale::Union{Nothing, Tuple{Symbol, Symbol}}
end

sname(s, comp) = Symbol("pb_", s, "_", comp)
rname(comp) = Symbol("R_", comp)

"All reactions of the model (rate keys are resolved per compartment, e.g. :kcleave_central)."
function reactions()
    rx = Rxn[]
    for comp in COMPS
        S(s) = sname(s, comp)
        R = rname(comp)
        kcl = Symbol("kcleave_", comp)
        for s in SPEC
            a = arms(s)
            for i in 1:2
                x, y = a[i], a[3 - i]
                if x == :c
                    push!(rx, Rxn([S(s)], [S(canon(:o, y))], :kopen, 1, nothing))
                    push!(rx, Rxn([S(s)], [S(canon(:m, y))], kcl, 1, nothing))
                elseif x == :o
                    push!(rx, Rxn([S(s)], [S(canon(:c, y))], :kclose, 1, nothing))
                    push!(rx, Rxn([S(s)], [S(canon(:m, y))], kcl, 1, nothing))
                    push!(rx, Rxn([R, S(s)], [S(canon(:oR, y))], isbound(s) ? :kon2 : :kon1, 1, nothing))
                elseif x == :m
                    push!(rx, Rxn([R, S(s)], [S(canon(:mR, y))], isbound(s) ? :kon2 : :kon1, 1, nothing))
                elseif x == :oR
                    push!(rx, Rxn([S(s)], [R, S(canon(:o, y))], :koff1, 1, nothing))
                    push!(rx, Rxn([S(s)], [S(canon(:mR, y))], kcl, 1, nothing))
                elseif x == :mR
                    push!(rx, Rxn([S(s)], [R, S(canon(:m, y))], :koff1, 1, nothing))
                end
            end
            if isbound(s)   # internalisation (central: loss only, as in the supplement matrix)
                if comp === :central
                    push!(rx, Rxn([S(s)], Symbol[], :kendo, 1, nothing))
                else
                    push!(rx, Rxn([S(s)], [Symbol("pb_", hasmR(s) ? :tmdd_m : :tmdd_o, "_", comp)], :kendo, 1, nothing))
                end
            elseif comp !== :tumor   # elimination of unbound forms
                push!(rx, Rxn([S(s)], [Symbol("pb_el_", comp)], :kel, 1, nothing))
            end
        end
        push!(rx, Rxn([R], Symbol[], :kendo, 1, nothing))
    end
    for s in UNBOUND   # exchange of unbound forms
        push!(rx, Rxn([sname(s, :central)], [sname(s, :peripheral)], :k12, 1, (:V1, :V2)))
        push!(rx, Rxn([sname(s, :peripheral)], [sname(s, :central)], :k21, 1, (:V2, :V1)))
        push!(rx, Rxn([sname(s, :central)], [sname(s, :tumor)], :k13, 1, (:V1, :V3)))
        push!(rx, Rxn([sname(s, :tumor)], [sname(s, :central)], :k31, 1, (:V3, :V1)))
    end
    return rx
end

"Ordered list of all state names (R, unbound, bound, internalised/eliminated pools) per compartment."
function state_names()
    out = Symbol[]
    for comp in COMPS
        push!(out, rname(comp))
        append!(out, [sname(s, comp) for s in SPEC])
        comp === :central && push!(out, :pb_el_central)
        if comp !== :central
            push!(out, Symbol("pb_tmdd_m_", comp), Symbol("pb_tmdd_o_", comp))
            comp === :peripheral && push!(out, :pb_el_peripheral)
        end
    end
    return out
end

# rate key of a Rxn resolved per compartment (kel, kendo, kon*, koff1, kopen, kclose are compartment independent)
const PARAM_KEYS = [:day, :kopen, :kclose, :kcleave_central, :kcleave_peripheral, :kcleave_tumor, :kel, :kendo, :kon1, :kon2, :koff1,
    :k12, :k21, :k13, :k31, :V1, :V2, :V3, :ksyn_peripheral, :ksyn_tumor]

# --- parameter sets ------------------------------------------------------------------------------
const DAY = 86400.0          # s per day: ODE time is in days, rate constants stay in 1/s as printed
const MW = 150_000.0         # g/mol, NOT STATED in the paper; 150 kDa IgG. Checked: 4.5 mg/kg*70 kg -> 808 nM in 2.6 L (Fig 4a, t = 0)
const KOPEN = 1.16e-2        # 1/s, supplement s003 parameter list entry 1 (kclose = Kmask * kopen; supplement 0.232 = 20 * 1.16e-2)

"Cynomolgus monkey parameters, Table 1 (rate constants 1/s, volumes L, synthesis nmol/s)."
function monkey_pars(; Kmask, kcleave, kopen = KOPEN)
    return Dict{Symbol, Float64}(
        :day => DAY, :kopen => kopen, :kclose => Kmask * kopen,           # Table 1: fold-masking Kmask; split of kopen/kclose: supplement
        :kcleave_central => kcleave, :kcleave_peripheral => kcleave, :kcleave_tumor => 0.0,  # Table 1: <3e-7; no tumour in monkey
        :kel => 1.0e-6,                                                   # Table 1
        :kendo => 1.0e-4,                                                 # Table 1
        :kon1 => 1.0e-3, :kon2 => 1.0e-3, :koff1 => 2.0e-3,               # Table 1
        :k12 => 1.1e-5, :k21 => 1.0e-5,                                   # Table 1
        :k13 => 0.0, :k31 => 0.0,                                         # no tumour in monkey
        :V1 => 0.1, :V2 => 0.1, :V3 => 1.0,                               # Table 1 (V3 unused)
        :ksyn_peripheral => 9.0e-5, :ksyn_tumor => 0.0)                   # Table 1 (nmol/s)
end

"Human projection parameters, Table 2. `kcleave_tumor_fold` = tumour kcleave relative to the systemic kcleave."
function human_pars(; Kmask, kcleave = 3.0e-7, kcleave_tumor_fold = 10.0, kopen = KOPEN, V2 = 2.6)
    return Dict{Symbol, Float64}(
        :day => DAY, :kopen => kopen, :kclose => Kmask * kopen,
        :kcleave_central => kcleave, :kcleave_peripheral => kcleave, :kcleave_tumor => kcleave * kcleave_tumor_fold,
        :kel => 6.3e-7,                                                   # Table 2
        :kendo => 1.0e-4,                                                 # Table 2
        :kon1 => 1.0e-3, :kon2 => 1.0e-3, :koff1 => 2.0e-3,               # Table 2
        :k12 => 4.8e-6, :k21 => 4.4e-6,                                   # Table 2
        :k13 => 1.9e-8, :k31 => 1.0e-5,                                   # Table 2
        :V1 => 2.6, :V2 => V2, :V3 => 0.01,                               # Table 2; V2 not printed (see fixes / notes)
        :ksyn_peripheral => 2.3e-3, :ksyn_tumor => 2.7e-4)                # Table 2 (nmol/s)
end

dose_nM(dose_mgkg, bw, V1) = dose_mgkg * bw * 1.0e6 / MW / V1   # mg/kg * kg -> mg -> nmol -> nM

# --- MTK system ----------------------------------------------------------------------------------
"""
    build_system()

Compile the 55-state ODE system once; parameters (mask strength, cleavage, species) are set via `ODEProblem`.
Observed quantities: `Ctot_*` (nM, all unbound and receptor-bound forms), `Cfree_*` (nM, unbound forms only) and
`Upt_peripheral/tumor` (nmol, cumulative internalised pools).
"""
function build_system()
    names = state_names()
    X = Dict{Symbol, Any}()
    for n in names
        X[n] = only(@variables $n(t) = 0.0)
    end
    P = Dict{Symbol, Any}()
    for k in PARAM_KEYS
        P[k] = only(@parameters $k = 1.0)
    end
    rhs = Dict{Symbol, Any}(n => Num(0) for n in names)
    for rx in reactions()
        rate = P[rx.rate]
        flux = rx.mult * rate * prod(X[s] for s in rx.from)
        for s in rx.from
            rhs[s] -= flux
        end
        gain = rx.vscale === nothing ? flux : flux * P[rx.vscale[1]] / P[rx.vscale[2]]
        for s in rx.to
            rhs[s] += gain
        end
    end
    # target synthesis (zero order): ksyn [nmol/s] / V [L] = nM/s
    rhs[:R_peripheral] += P[:ksyn_peripheral] / P[:V2]
    rhs[:R_tumor] += P[:ksyn_tumor] / P[:V3]
    eqs = Equation[D(X[n]) ~ P[:day] * rhs[n] for n in names]
    tot(comp) = sum(X[sname(s, comp)] for s in SPEC)
    free(comp) = sum(X[sname(s, comp)] for s in UNBOUND)
    @variables Ctot_central(t) Ctot_peripheral(t) Ctot_tumor(t) Upt_peripheral(t) Upt_tumor(t)
    @variables Cfree_central(t) Cfree_peripheral(t) Cfree_tumor(t)
    append!(eqs, [Cfree_central ~ free(:central), Cfree_peripheral ~ free(:peripheral), Cfree_tumor ~ free(:tumor),Ctot_central ~ tot(:central), Ctot_peripheral ~ tot(:peripheral), Ctot_tumor ~ tot(:tumor),
        Upt_peripheral ~ (X[:pb_tmdd_m_peripheral] + X[:pb_tmdd_o_peripheral]) * P[:V2],
        Upt_tumor ~ (X[:pb_tmdd_m_tumor] + X[:pb_tmdd_o_tumor]) * P[:V3]])
    @named stroh = System(eqs, t)
    return mtkcompile(stroh)
end

"Initial-condition and parameter pairs for `ODEProblem(sys, pairs, tspan)`: baseline receptor at steady state."
function problem_pairs(sys, pars::Dict{Symbol, Float64})
    pr = [getproperty(sys, k) => v for (k, v) in pars]
    ic = [sys.R_peripheral => pars[:ksyn_peripheral] / pars[:kendo] / pars[:V2],     # steady state ksyn/(kendo*V)
          sys.R_tumor => pars[:ksyn_tumor] / pars[:kendo] / pars[:V3]]
    return vcat(pr, ic)
end
