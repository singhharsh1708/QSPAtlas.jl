# Validation

Two checks back the atlas: the population simulation layer was compared with an
established open-source tool, and every model was implemented a second time from its paper
alone, with separate code, and compared with the atlas implementation. Both checks were
made within this project; none has yet been made by outside people.

## Simulation layer against mrgsolve

The dosing and population machinery was run against [mrgsolve](https://mrgsolve.org)
2.0.1 (R) on identical models, doses and time grids, with tight solver tolerances on both
sides (2026-11-04). The scripts are in
[`validation/mrgsolve/`](https://github.com/singhharsh1708/QSPAtlas.jl/tree/main/validation/mrgsolve).

| case | largest relative difference | result |
|---|---|---|
| IV bolus, repeated doses (`ii`, `addl`) | 8.7e-10 | agree |
| infusions every 12 h, two overlapping infusion streams, infusion longer than `ii` | 5.4e-10 | agree (rates add) |
| oral dose with bioavailability 0.6 and lag 0.75 h; lag longer than `ii` | 7.4e-9 | agree (mrgsolve fails on lag longer than `ii` with `addl`; explicit dose records agree) |
| bioavailability on a fixed-rate infusion (scales duration) and on a fixed-duration infusion (scales rate) | 5.7e-10 | agree |
| steady state: bolus, oral with lag, infusion shorter and longer than `ii`, turnover model, steady state after earlier doses | 7.2e-9 | agree |
| steady state with lag longer than `ii` | 4e-9 against 200 explicit doses | QSPAtlas agrees with the explicit doses; mrgsolve fails |
| reset (EVID 3) and reset with dose (EVID 4) during an infusion | 3.3e-10 | agree |
| new infusion after a reset | | QSPAtlas matches the closed form (17.470 at 12 h); mrgsolve keeps the old infusion's stop event and ends the new infusion 6 h early |
| doses at the same time; a bolus at the moment one infusion stops and another starts | 5.5e-10 | agree |
| population: one-compartment oral, correlated log-normal CL and V, proportional error, 5,000 and 20,000 subjects | within Monte Carlo error | agree; the zero-random-effect prediction equals the deterministic model in both tools |

### Dataset cases

NONMEM-style dataset files (`validation/mrgsolve/datasets/`) were read by both tools,
`read_nmdata` in QSPAtlas and `mrgsim(data = ...)` in mrgsolve, with mrgsolve's default record
handling (data set order kept, `recsort = 1`, `nocb = TRUE`) on a two-compartment model with
depot and a turnover response, clearance and central volume scaled by a covariate WT
(2026-11-05). One row per observation record:

| case | largest relative difference | result |
|---|---|---|
| record order: observations before and after dose records at the same time (bolus, oral, infusion start), observation at an `ADDL` dose time | 4.7e-9 | agree |
| time-varying covariate, NOCB | 1.4e-9 | agree |
| time-varying covariate, LOCF | 1.4e-9 (states) | states agree; the output at a record where the covariate changes differs (see below) |
| EVID 2 records carrying covariate changes | 1.4e-10 | agree; no output rows for EVID 2 |
| SS=2: 10 every 24 h at 8 h plus 15 at 20 h; 10 every 12 h plus 5 every 24 h; oral SS=2 on an infusion SS=1 | 6.1e-9 (PK) | PK agrees; the turnover response differs (see below) |
| steady-state constant infusion (AMT 0, RATE 30, SS 1), alone, continued by an explicit infusion, and with a modelled rate (RATE -1) | 1.7e-10 | agree; in both tools the infusion ends at the record, as the NONMEM Users Guide describes |
| modelled rate (RATE -1) and duration (RATE -2) with bioavailability 0.5, also with SS=1 | 2.6e-10 | agree |
| reset (EVID 3) and reset with dose (EVID 4) with observations before and after them | 6.1e-9 | agree |
| SS=1 record while `ADDL` doses of an earlier record are still to come (bolus doses; infusions, one running at the SS=1 time) | 6.9e-10 on the 14 records the cancellation does not reach (before the first cancelled dose, and a control subject without SS); up to 0.93 after it | differ after the SS=1 record: mrgsolve keeps the pending `ADDL` doses (see below) |

The same semantics are tested against closed forms in `test/dataset.jl` (SS=2 against
superposed and explicit dose histories, the constant infusion against `Css = R / CL`).

### Dose control with random effects and between-occasion variability

A one-compartment model with first-order absorption, bioavailability
`F = 0.6 exp(ηF) WT / 70` and lag `ALAG = 0.75 exp(ηLAG)` on the depot (computed by the model's
`dosecontrol`), and clearance `CL = 2 exp(κ)` with κ of the occasion (OCC column, NOCB; mrgsolve
selects κ by OCC in `$MAIN`, the NONMEM indicator-variable way), 2026-11-05. In the exact case
every subject has SUPPLIED η and κ values (data set columns in mrgsolve, `η` and `κ` arguments
in QSPAtlas), so individual curves are compared record by record; in the distribution case
both tools draw their own random effects.

| case | largest relative difference | result |
|---|---|---|
| fixed η and κ: two occasions, doses at the start of each | 6.2e-11 | agree |
| fixed η and κ: WT changes on the second dose record (that dose uses WT 95) | 1.1e-10 | agree |
| fixed η and κ: `addl` doses, and `addl` doses running into occasion 2 and a new WT (they keep their record's F and ALAG in both tools) | 3.8e-11 | agree |
| fixed η and κ: occasion labels 1, 2, 1 (the first κ comes back) | 9.0e-11 | agree |
| sampled η (F, ALAG) and κ (CL, two occasions), 5,000 subjects: 5th, 50th, 95th percentiles at 14 times | largest difference 1.77 bootstrap standard errors (none above 2) | agree within Monte Carlo error |

The same features are tested against closed forms in `test/iov.jl` (log-normal AUC under
log-normal F, the lag distribution, covariance and independence of the κ draws, exact switching
at occasion boundaries).

For the population case both tools were also compared with a near-exact reference (closed
form, 4 million draws). QSPAtlas sampling with one million draws matches the reference to
1.7e-3 (largest z-score 1.96).

Differences in convention, documented in [Population simulation](simulation.md):

- **Where steady state starts.** QSPAtlas iterates dosing cycles from the model's initial
  conditions, mrgsolve from the state at the steady-state record. QSPAtlas and mrgsolve agree
  on every steady-state case above; they can differ for models with several periodic
  solutions or slow, history-dependent states. How NONMEM initialises its steady-state
  computation depends on the routine and was not compared with NONMEM output here.
- **Observations at a dose time** see the state after the dose (mrgsolve `recsort = 3`)
  unless tagged `:pre`; datasets set the tag from the record order, and `ADDL` doses come
  before an observation at the same time (mrgsolve's default; not checked against NONMEM).
- **LOCF covariates.** With `nocb = FALSE`, mrgsolve evaluates the output of a record with the
  previous record's covariates (the record's own values start after it); QSPAtlas uses the
  record's values at the record. The states agree; an output that depends on the covariate
  (here the concentration, through the volume) differs at the change records. NOCB is the
  convention the NONMEM documentation describes; under NOCB QSPAtlas and mrgsolve agree.
- **SS=2 with an endogenous baseline.** mrgsolve adds the whole steady state of the SS=2 dose,
  including the turnover response's own baseline (the response jumps from 100 to 200 with no
  drug effect when IMAX = 0). QSPAtlas adds the periodic state minus the dose-free periodic
  state, so the baseline is not counted twice. For PK compartments starting at zero the two
  are the same and agree to 6e-9.
- **Steady-state constant infusion.** The infusion ends at the record (NONMEM Users Guide);
  mrgsolve 2.0.1 behaves the same, although its user guide describes the infusion as never-ending.
- **SS=1 after a record with `ADDL` doses.** QSPAtlas cancels the `ADDL` doses of earlier
  records that are still to come at the SS=1 time, as the NONMEM help on SS describes
  ("infusions in progress or pending additional doses are cancelled"); not compared with
  NONMEM output. mrgsolve 2.0.1 stops the running infusion but keeps the later `ADDL` doses:
  its output agrees to 7.0e-10 with QSPAtlas run with those doses written as dose records of
  their own (case `d8_ss_addl`, `ss_addl_qspatlas.jl` and `ss_addl_mrgsolve.R`). QSPAtlas on
  the dataset and QSPAtlas with the cancelled doses left out by hand agree exactly.
- **A reset** cancels lagged doses that have not yet been given and keeps later `addl` doses.
- **A lagged dose entering exactly at an observation time.** QSPAtlas reports the state after
  the dose (as for any dose at an observation time); mrgsolve 2.0.1 reports the state before it.
  Found while setting up the IOV cases (a lag of exactly 0.75 h on a 0.75 h observation); the
  cases above use lags that do not coincide with observation times.

## Third implementation from the authors' own code

For four models the authors published their NONMEM or Monolix code. That code was
translated to [rxode2](https://nlmixr2.github.io/rxode2/) (R) with nonmem2rx and
monolix2rx (or by hand where the code is printed in a PDF), run with the parameter values
and inputs the atlas uses (the code often carries initial estimates, so this checks
equations, structure, dosing and units, not the values), and compared point by point with
the atlas. Scripts:
[`validation/authors_code/`](https://github.com/singhharsh1708/QSPAtlas.jl/tree/main/validation/authors_code).

| model | route | cases | largest relative difference |
|---|---|---|---|
| Ouerdani 2015 | nonmem2rx (both control streams) | Figure 4 mouse and patient; mouse AUC 220 and 657 | 5.2e-10 |
| Hutchinson 2016 | monolix2rx and a hand translation (identical) | Figures 2a, 2c, 2e and 6a | 1.7e-10 (2.5e-10 in the window cases with a tighter atlas solve) |
| Lindauer 2017 | nonmem2rx, unedited | mouse PK, receptor occupancy and tumour volume at 0, 0.1, 1.4 and 5 mg/kg | 1.2e-8 |
| Stein 2019 | by hand (code printed in the supplement PDF) | typical curve and three comedication orders | 2.2e-12 |

Differences between the code and the paper that this check exposed (none affects an atlas
target): in Ouerdani's code a control run (DOSE = 0) switches all drug terms off, while the
paper's Eq. 2 with beta_alpha = 0 keeps alpha = alpha0 for any exposure; in Hutchinson's code
the normalisation factor is N = Nmax x dose + 1, not Eq. 5's Nmax; Hutchinson's Figure 6
chemotherapy is not in the code. Stroh and Fornari published no runnable code (Stroh's
supplement has the reaction matrices, checked entry by entry; Fornari's code is truncated).

## Second implementations of the six models

Each model was implemented a second time from its paper and supplement alone, without
access to the atlas code, in plain OrdinaryDiffEq without ModelingToolkit, within this
project. The second implementation made its own choices for every unstated input. It
compared simulations only: nothing was re-digitised, so digitising errors are not covered
by this check. Where both implementations chose an unstated input (a dose group, a starting
value) by matching the figure, their agreement is not independent evidence for that choice.
Results:

| model | agreement with the atlas implementation |
|---|---|
| Ouerdani 2015 | Figure 4 curves agree to 0.003 mm³ (mice) and 0.06 mm (patient K); the patient P curves differ by 0.2 mm at the start only, because the two implementations read the starting size from the figure as 66.8 and 67 mm (66.8 vs 67 mm). The same dose group (100 mg/kg) was chosen independently. |
| Hutchinson 2016 | All 22 Figure 2 curves agree to within 2e-5 of each curve's range. The rebuild found a reading of Figure 6 (resistance clock from the start of chemotherapy, chemotherapy scale 0.01) that brings panels b to d close; without independent evidence for it, it is reported as a sensitivity run and Figure 6 stays C. |
| Lindauer 2017 | Identical base values (-1.302 and -0.401 with the Table 1 KdegPD1) and the same six bars that do not match the figure, with the same unstated choices made independently. |
| Fornari 2019 | Same closure fixes. The rebuild found that reading the rat PK volumes per kg reproduces Figure 3 (see the model page). |
| Stroh 2019 | The reaction network built from the text matches the supplement matrices entry for entry (234 and 84 entries). |
| Stein 2019 | Same reading of the variability (standard deviations) and the same result: early phase inside the published bands, late phase below them in every variant tried. |

No transcription error was found in any atlas model. One second implementation (Fornari)
changed the reading of a units question that decides a figure grade; see that model's page.
