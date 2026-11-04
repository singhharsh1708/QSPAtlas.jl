# stroh_2019 research notes

Stroh M, Sagert J, Burke JM, Apgar JF, Lin L, Miliard BL, Kavanaugh WM. Quantitative Systems Pharmacology Model of a Masked, Tumor-Activated Antibody. CPT Pharmacometrics Syst Pharmacol 2019;8:676-684. DOI 10.1002/psp4.12448, PMC6765697.

## Sources

- Main text Tables 1 (monkey) and 2 (human); Fig 2a/b schematic; Fig 3 (monkey), Fig 4 (human), Fig 5 (multi-dose).
- Supplement s003 = KroneckerBio matrix code: 55 states, dx/dt = k + A1 x + A2 (x kron x). The printed "k" is NOT a zero-order vector: it is a 32-entry parameter list (kopen 1.16e-2, kclose 0.232, kcleave 1e-6 / 1e-6 / 1e-5, kel 6.3e-7, kon 1e-3, koff 2e-3, kendo 1e-4, 2*kopen, 2*kclose, 2*kcleave, ksynR 2.34e-3 and 2.8221e-5, k21 4.4e-6, k12 4.8e-6, k31 4.4e-6, k13 2.5e-9). The true zero-order vector (receptor synthesis) is not printed. The snapshot is one human-like parameter set (Kmask 20, kcleave_central 1e-6) that does NOT equal Table 2.
- State layout (verified from the index pattern, N = 55, A2 column = (i-1)*55 + j): 1 = R_central, 2-17 central (6 unbound c_c, c_o, o_o, c_m, o_m, m_m; 9 receptor-bound; el), 18 = R_periph, 19-36 periph (adds tmdd_m, tmdd_o, el), 37 = R_tumour, 38-55 tumour (no el state used).
- s001 = Fig S1 (likelihood scans), s002 = Fig S2 (cleaved flux): not needed.

## Targets: pre-declared (before any simulation was compared)

Grade-bearing (all model lines, kind: curve):

- Figure 3a (3 mg/kg: mAb0,0, Pb-Tx M1,S1, M2,S2, M1,0), 3b (5 mg/kg: mAb0,0, M2,S1, M1,S2), 3c (10 mg/kg: mAb0,0, M2,S1, M1,S2): 10 curves, log-y plasma total drug [nM]. (Fig 3 lines start at day about 1-10 where the digitiser could separate them from circles and other curves; the span actually traced is in each CSV.)
- Figure 4a-e: 25 curves (mAb0,0, 10x, 100x, 1000x, 10000x mask) for plasma, periphery, tumour total drug (log y) and cumulative receptor-mediated uptake in periphery / tumour (linear y). Overlapping blue / orange curves are only graded over the span where they were traceable.

Not graded: observed monkey points (hollow circles), Figure 5 (multiple dose; optional, not digitised).

Tolerance: default (rms 0.015, p95 0.03 of the axis box). No reason to deviate was recorded.

## Modelling decisions made before comparing (written up as assumptions / fixes)

1. Dose: MW 150 kDa is not stated in the paper. Check: 4.5 mg/kg * 70 kg / 150 kDa / 2.6 L = 808 nM, matching Fig 4a at t = 0 (about 800 nM).
2. Breathing absolute rates are not in the main text (only Kmask = kclose/kopen). Primary reading: kopen = 1.16e-2 1/s (supplement parameter 1), kclose = Kmask * kopen (supplement: 0.232 = 20 * 1.16e-2). The alternative reading (kclose = 0.232 fixed, kopen = 0.232/Kmask) is run as a sensitivity only.
3. The dose enters as intact c_c (both masks closed) for Pb-Tx, as m_m for the parental mAb0,0.
4. kcleave = 3e-7 1/s for S1 and S2 (upper bound, Table 1), 0 for Pb-Tx M1,0 (no substrate); mAb0,0 is already fully cleaved. Cleavage per masked site: c_c -> c_m at 2 kcleave etc. (checked against s003).
5. Compartments: transfer rates are amount based (the k13 = Q p/(p + V1/V3) derivation gives C3/C1 = p), implemented with volume ratios.
6. Human: V2 is not in Table 2; inferred 2.6 L (ksynR_periph 2.3e-3 nmol/s = 9 nM monkey concentration * 2.6 L * kendo). kcleave in tumour = 10x central (text: "10-fold increased kcleave in the tumor" is the scenario shown in Fig 4d, e; the supplement snapshot also has 10x). Tumour ksynR = 2.7e-4 nmol/s (Table 2; supplement parameter 2.8221e-5 conflicts: recorded).
7. Plotted "total drug" = sum of all unbound and receptor-bound species in the compartment (excluding internalised / eliminated pools); "receptor-mediated uptake" = cumulative tmdd_m + tmdd_o pools (nmol).
8. The central receptor R_central is in the supplement state list but has no stated expression: initial 0, no synthesis.

## Freeze

Targets digitised and written to `reference/*.csv` and `reference/digitise.toml`: 2026-11-04 09:48 (local), before any simulation of this model was run. `protocol.jl` TOLERANCE set to the defaults at that time.

## Outcome (2026-11-04)

- Matrix structure check passed: 234/234 entries of A1 and 84/84 of A2 identical to the transcribed supplement matrices (`supplement_check/check_matrices.jl`).
- Grade C. As published, 14/35 curves pass (Fig 3: 7/10; Fig 4: 7/25); with the `unbound_readout` from_target fix, 18/35.
- Observed after comparing (not before): periphery / tumour "Total Drug" of the mAb curves only matches as unbound forms (evidence: 1e-2 nM tumour mAb << R_T); strongly masked curves clear faster in the model than in the figure; the sensitivity scan (`sensitivity.jl`, table in `report.md`) shows that kcleave of about 1e-7 would give 24/25 Fig 4 passes. This is a parameter choice by figure match, so it is reported only as a sensitivity, not as a fix and not in the grade.
- The breathing split (kopen vs kclose) was pre-declared as kopen fixed; the alternative runs show it is immaterial.

## Review (2026-11-04, one round)

MW 150 kDa and human V2 2.6 L were made explicit :missing Fix entries (both runs); the kcleave wording was softened; the Fig 3 traces keep t >= 1 d (min_x rule; contaminated marker points, including one at negative time, were removed; no pass / fail change); the readout stays :from_target. Final grade after review: C.
