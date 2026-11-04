# Validation

Two independent checks back the atlas: the population simulation layer was compared with
an established open-source tool, and every model was rebuilt a second time from its paper
alone, with separate code, and compared with the atlas implementation.

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

For the population case both tools were also compared with a near-exact reference (closed
form, 4 million draws). QSPAtlas sampling with one million draws matches the reference to
1.7e-3 (largest z-score 1.96).

Differences in convention, documented in [Population simulation](simulation.md):

- **Where steady state starts.** QSPAtlas iterates dosing cycles from the model's initial
  conditions, mrgsolve from the state at the steady-state record, NONMEM from zero. All
  three agree when the iteration converges; they differ only for models with several
  periodic solutions or slow, history-dependent states.
- **Observations at a dose time** see the state after the dose (mrgsolve `recsort = 3`).
- **A reset** cancels lagged doses that have not yet been given and keeps later `addl` doses.

## Independent rebuilds of the six models

Each model was rebuilt from its paper and supplement alone, without access to the atlas
code, in plain OrdinaryDiffEq without ModelingToolkit. The rebuild made its own choices for
every unstated input. Results:

| model | agreement with the atlas implementation |
|---|---|
| Ouerdani 2015 | Figure 4 curves agree to 0.003 mm³ (mice) and 0.06 mm (patients); the only difference is the starting size read from the figure (66.8 vs 67 mm). The same dose group (100 mg/kg) was chosen independently. |
| Hutchinson 2016 | All 22 Figure 2 curves agree to within 2e-5 of each curve's range. The rebuild found a reading of Figure 6 that reproduces panels b to d (see the model page). |
| Lindauer 2017 | Identical base values (-1.302 and -0.401 with the Table 1 KdegPD1) and the same six bars that do not match the figure, with the same unstated choices made independently. |
| Fornari 2019 | Same closure fixes. The rebuild found that reading the rat PK volumes per kg reproduces Figure 3 (see the model page). |
| Stroh 2019 | The reaction network built from the text matches the supplement matrices entry for entry (234 and 84 entries). |
| Stein 2019 | Same reading of the variability (standard deviations) and the same result: early phase inside the published bands, late phase below them in every variant tried. |

No transcription error was found in any atlas model.
