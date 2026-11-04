# fornari_2019 working notes

Fornari C, O'Connor LO, Pin C, Smith A, Yates JWT, Cheung SYA, Jodrell DI, Mettetal JT, Collins TA. Quantifying Drug-Induced
Bone Marrow Toxicity Using a Novel Haematopoiesis Systems Pharmacology Model. CPT Pharmacometrics Syst Pharmacol 2019;8:858-868.
DOI 10.1002/psp4.12459, PMC6875710. (Author list checked against the supplement front page: the paper lists
Fornari, Oplustil O'Connor, Pin, Smith, Yates, Cheung, Jodrell, Mettetal, Collins; the supplement writes "Lenka Oplustil O'Connor".)

## What is where
- Main text Eqs 1-13 (MathML in the article XML, converted to text). Tables 1, 2 (rat + human parameters).
- Supplement s012 (text S1-S12): PK model S5 (Eq S1), steady-state algebra S6/S7 (Eq S3-S7), human scaling S9-S12.
- s013 "Supplementary Model Code": PDF has 4 pages; it contains the rat parameter script (all rat values incl. alphaWBC = 0.5 and the
  parameter mapping of the ODE function) but the page ends right after `k12_Carbo = QC/Vcen_Carbo;`. The ODE right-hand sides, the dosing
  call and the solver call are NOT in the file (truncated). Checked by rendering page 4 as an image.
- s008 Table S1 (rat study designs), s009 Table S2 (IC50), s010 Table S3 (rat PK), s011 Table S4 (residual error).
- s001 Fig S1 (30 mg/kg one cycle), s003 Fig S3 (rat free-Pt PK + Damage, vector PDF, no embedded image).

## Equations as printed (and problems)
Eq 1  dMPP/dt = kstem fdbk_stem + kprol_MPP (1 - Slope_MPP Damage) MPP - ktr_CMP fdbk_mat_MPP - ktr_MEP fdbk_mat_MPP
      (printed `fdbk_mat_MPP` does not exist in Eq 12; maturation terms lack the factor MPP. Eq 4 gain terms and S6 steady state
      show the intended form: - ktr_CMP fdbk_mat_CMP MPP - ktr_MEP fdbk_mat_MEP MPP.)
Eq 2  drug effect = Slope x Damage x kprol x Cells (kills the proliferating flux; (1 - Slope Damage) can become negative)
Eq 3  dDamage/dt = kdam Cp_freePt - krep Damage
Eq 4  dCMP/dt = kprol_CMP (1 - Slope_CMP Damage) CMP + ktr_CMP fdbk_mat_CMP MPP - ktr_Neut CMP - ktr_Mono CMP
      dMEP/dt = kprol_MEP (1 - Slope_MEP Damage) MEP + ktr_MEP fdbk_mat_MEP MPP - ktr_Ret MEP - ktr_Plt MEP   (printed left side "dCMP/dt" twice: typo)
Eq 5  neut/mono 3 transit compartments: dT1 = ktr CMP - a T1; dT2 = a T1 - a T2; dT3 = a T2 - a T3
Eq 6  a_Neut = 3/MTT_Neut, a_Mono = 3/MTT_Mono
Eq 7  Ret/Plt: dT1 = kprol fdbk_prol T1 - delta Slope_MEP kprol Damage T1 + ktr MEP - a T1; dT2 likewise + a T1 - a T2; dT3 = a T2 - a T3
      (printed platelet T1 uses a_Ret: typo, should be a_Plt; "T1_Ptt" typo)
Eq 8  a_Ret = 3/MTT_Ret, a_Plt = 3/MTT_Plt
Eq 9  dNeut = a T3 - kcirc Neut; dMono; dRet; dPlt. dRBC is NOT printed (Table 1 has kcirc_RBC; S4 has kcirc_Ret = kcirc_RBC RBC0/Ret0).
Eq 10 Slope_H = Slope_R (IC50_H / IC50_R)^-1
Eq 11 fdbk_stem = (1/2 CMP0/CMP + 1/2 MEP0/MEP)^gamma_stem
Eq 12 fdbk_mat_CMP = (1/2 Neut0/Neut + 1/2 Mono0/Mono)^gamma_matCMP; fdbk_mat_MEP = (RBC0/RBC)^gamma_matMEP
Eq 13 fdbk_prol_Ret = (1/2 Ret0/Ret + 1/2 RBC0/RBC)^gamma_prolTrans; fdbk_prol_Plt = (Plt0/Plt)^gamma_prolTrans
      (the 1/2 weights are alphaWBC = 0.5 in the supplementary code parameter list)
PK (S5, Eq S1): dCen = -CL/Vcen Cen - Q/Vcen Cen + Q/Vper Per; dPer = Q/Vcen Cen - Q/Vper Per; Cp = Cen/Vcen. Table S3: QC 0.06, Vcen 0.26, Vper 0.30, CL 0.42 (units printed "1/h" for QC, CL; code says liter/hour).
Steady state (S3, S4, S7): all rate constants (ktr*, kprol*, kstem, a, kcirc_Ret, T1/T2 baselines) are DERIVED from baselines, MTTs,
kcirc and lambda1..lambda5 (Table 1). They are derived here directly from the steady-state equations (same results as S4 where legible).

## Missing / ambiguous inputs found
1. Rat dose -> amount conversion (mg/kg to umol in Cen) not stated anywhere. Body weight only "minimum 320 g" (S1); molecular weight not given.
   Fig S3 (PK model fit, 30 mg/kg) shows Cp(0) ~ 450 umol/L and a Damage peak ~4.77 => ~117 umol per 30 mg/kg (3.9 umol per mg/kg), which is
   4.5x what 0.32 kg x 371.25 g/mol gives (25.9 umol). Plan: as-published = 0.32 kg and MW 371.25 (the only stated/standard numbers);
   fix `dose_scale` = value inferred from Fig S3 (PK-only evidence, independent of Fig 3). (Revised after review, see the last section.)
2. ODE structural typos (above) - resolved by steady-state consistency (S6) - Fix kind :typo.
3. RBC equation absent - Fix :missing from S4 (kcirc_Ret definition).
4. Dose timing: Table S1 'dosing day 1'; Fig S3 x axis (hours) shows dose at t = 24 h; Fig 3 shows 8 doses (day 1,15,...,99; Table S1 lists
   1,15,29,43,57,71,85 for study E, garbled in the PDF; the text says eight cycles). Use dose at t = 24 h + 336 h*(i-1), i = 1..8.
5. Human baselines are RANGES in Table 1: Neut 2-7, Mono 0.2-10, Ret 40-115, RBC 4100-5900 (10^9/L). Only Plt 345, MPP 1.3, CMP 20.9,
   MEP 15 are single values. Fig 4 baseline neutrophils ~4.7 (red/blue lines at t = 0).
6. Human PK is not in the paper or supplement: S11 only says "2-compartment model of Zandvliet et al. (BJCP 2008;66:485, PMC2561111)"
   and Calvert dose = AUC (GFR+25) with AUC 5 mg min/mL, GFR 125 mL/min => 750 mg. Zandvliet Table 2 (read from PMC2561111):
   CL = 0.76 CrCl + 1.5 L/h, V1 15.5 L, Q 3.46 L/h, V2 9.86 L, ultrafiltrable Pt (umol/L), 30 min infusion in that study.
   Used external parameters as Fix :missing; unclear if GFR 125 mL/min = CrCl in the CL formula (assumed).
7. Eq 10 vs Table 2: Slope_R x IC50_R/IC50_H with Table S2 (CMP 2.83/4.04, MEP 0.86/1.56) gives CMP 1.03, MEP 1.21; Table 2 prints 0.57 and 0.66
   (and 0.79 for MPP, no MPP IC50 in S2). Printed Table 2 values are 0.388x (CMP) and 0.301x (MEP) of the rat slopes.
8. Fig 4 simulation details unspecified: individual variability (1000 virtual patients, no IIV in the paper), whether red line is typical value or mean,
   infusion duration, start time of dosing.
9. Fig 3 red line: "average model predictions" of 1000 simulations with the VPC; no IIV model given => typical-value simulation used.

## TARGETS, pre-declared (written before any simulation; digitisation done 2026-11-04 09:49 IST)
Graded (kind: curve, red mean-prediction lines), default tolerance (rms 0.015, p95 0.03 of the axis box) for all, no per-curve overrides:
 fig3_MPP, fig3_CMP, fig3_MEP (bone marrow), fig3_Neut, fig3_Plt (blood); fig4a_neut, fig4a_plt (rat slopes), fig4b_neut, fig4b_plt (IC50-scaled).
Not targeted (reported only): Fig 3 monocytes, reticulocytes, RBC (the lines are digitisable but not graded, to keep the target set manageable);
 Fig 3 observation points; Fig 4 blue Schmitt/Friberg lines and percentile bands (not the authors' model); Fig 4c bars; Supp figs S1, S2, S4-S7.
Calibration: see reference/digitise.toml. Log y axes calibrated from tick-label glyph centroids (-0.5 px) so +-0.5-1 px calibration uncertainty (0.4-0.8% of the box) is inside
the default tolerance but eats part of it; no loosening declared. Fig 3 Neut: only 54 of 136 columns traced (line is faint under data markers, early days missing);
Fig 4 red lines are partly hidden under the dark-blue Schmitt mean line (first ~3 days missing).
Verification overlay (local only, it contains the paper's figures): figures/digitise_check.png (Fig 3 above, Fig 4 below) (points sit on the red lines; extra ticks: Fig 3 x label 60, MEP 10^1; Fig 4 day-20 tick, threshold lines 1/0.5 and 50/25).
Human-parameter rule (declared before comparing): as-published run uses the midpoint of each printed range (Neut0 4.5, Mono0 5.1, Ret0 77.5, RBC0 5000), the Zandvliet
PK with 750 mg over 30 min, typical value, dose at t = 0, Table 2 human slopes as printed (b) and rat slopes (a). Sensitivity runs (not graded) afterwards.

## Dose scale evidence (Fig S3, PK-only; decided 2026-11-04 before any Fig 3/4 simulation)
`dose_scale_from_S3.jl`: with Table S3 PK parameters the digitised blue free-Pt line (clean
terminal phase t >= 27 h, 5.2 ... 0.12 uM) is reproduced within 1-3 % by 117 umol in Cen for 30 mg/kg (3.90 umol per mg/kg; terminal rate 0.1725/h also matches
the digitised 0.177/h, confirming the PK parameters). The orange Damage curve (peak 4.77) needs 132 umol (4.41 per mg/kg) with kdam = krep = 0.017, i.e. it is ~13 % higher
than the printed Damage ODE gives for the printed PK, and decays faster (0.0154/h vs 0.0127/h). Unexplained (issue to report).
DECISION (before comparing with Fig 3): fix dose_scale uses the free-Pt value 3.90 umol per mg/kg (PK curve is a direct fit of the PK model to concentration data; the Damage curve is a derived
output). The 4.41 variant and the physical 0.32 kg / 371.25 g/mol (0.862) variant are run as sensitivity only and do not change the grade.
Published run: cannot run as printed (Eq 1 refers to an undefined fdbk_mat_MPP and lacks the MPP factor, Eq 9 has no RBC equation); implemented as an
error in the as-published run, so grade is from the fixed run (B if all targets pass, else C).

## Results (2026-11-04, after freeze; first run, with dose_scale taken from Fig S3)
Grade C. With fixes: MPP 0.0105, MEP 0.0099, Plt 0.0075 pass; CMP 0.0316, Neut 0.0355 fail; all four Fig 4 lines fail (0.040-0.069). As published: error (not closed).
A bug found and fixed on the way: the dosing callback was passed both to ODEProblem and to solve (doses applied twice, Damage peak doubled); the first run's numbers were discarded.
Sensitivity runs (`sensitivity.jl`) are summarised in `report.md`, which also has the issues list.

## Review (2026-11-04)
The review graded the model D, because dose_scale had been set from Fig S3 (a value fitted to a figure). Done: dose_scale = 0.320 kg (S1 minimum weight) x 1000 / 371.25; S3 value 3.9 (~1.45 kg) diagnostic only. human_schedule fix (Fig 4 caption: two 21-day cycles) added. Rescored: Fig 3 all five fail (rms 0.047-0.134), Fig 4 passes only 4b neutrophils (0.0099); grade C. Report corrected.
