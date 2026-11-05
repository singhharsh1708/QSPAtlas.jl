# Cross-validation against mrgsolve

Scripts that run the same models, doses and time grids in QSPAtlas and in
[mrgsolve](https://mrgsolve.org) (R) and compare the results. Summary and verdicts:
[validation page](https://singhharsh1708.github.io/QSPAtlas.jl/dev/validation/).

| script | what it does |
|---|---|
| `cases_qspatlas.jl`, `cases_mrgsolve.R` | deterministic cases 1 to 7: repeated bolus doses, overlapping infusions, bioavailability and lag, F on fixed-rate and fixed-duration infusions, steady state, resets, simultaneous events |
| `compare.R` | per-case largest absolute and relative differences |
| `pop_qspatlas.jl`, `pop_mrgsolve.R`, `pop_compare.R` | case 8: population quantiles (5,000 and 20,000 subjects) with bootstrap intervals of the differences |
| `pop_reference.R` | near-exact closed-form reference for the population quantiles |
| `pop_rng_check.jl`, `pop_rng_check.R` | sampling checks of the random effects in both tools |
| `datasets/d1_*.csv` to `d7_*.csv`, `dataset_cases_qspatlas.jl`, `dataset_cases_mrgsolve.R`, `dataset_compare.R` | dataset cases d1 to d7: the same NONMEM-style files read by both tools (record order, NOCB and LOCF covariates, EVID 2, SS=2, steady-state constant infusion, RATE -1 / -2, EVID 3 / 4), compared record by record |
| `datasets/d8_ss_addl.csv`, `ss_addl_qspatlas.jl`, `ss_addl_mrgsolve.R` | dataset case d8: an SS=1 record while `ADDL` doses of an earlier record are still to come. QSPAtlas cancels them (NONMEM help, SS) and is also checked against the same subjects with the cancelled doses left out by hand; mrgsolve 2.0.1 keeps them, so the tools differ after the SS=1 record. `dataset_compare.R` picks up the outputs |
| `datasets/iov_exact.csv`, `iov_cases_qspatlas.jl`, `iov_cases_mrgsolve.R`, `iov_compare.R` | between-subject variability on F and ALAG (with a covariate on F) and between-occasion variability on CL (OCC column): supplied η and κ per subject compared record by record, and a 5,000-subject distribution check with sampled random effects |

Run from the repository root:

```bash
julia --startup-file=no --threads=4 --project=. validation/mrgsolve/cases_qspatlas.jl
julia --startup-file=no --threads=4 --project=. validation/mrgsolve/pop_qspatlas.jl
Rscript validation/mrgsolve/cases_mrgsolve.R && Rscript validation/mrgsolve/compare.R
Rscript validation/mrgsolve/pop_mrgsolve.R && Rscript validation/mrgsolve/pop_compare.R
julia --startup-file=no --threads=4 --project=. validation/mrgsolve/dataset_cases_qspatlas.jl
Rscript validation/mrgsolve/dataset_cases_mrgsolve.R && Rscript validation/mrgsolve/dataset_compare.R
julia --startup-file=no --threads=4 --project=. validation/mrgsolve/ss_addl_qspatlas.jl
Rscript validation/mrgsolve/ss_addl_mrgsolve.R && Rscript validation/mrgsolve/dataset_compare.R
julia --startup-file=no --threads=4 --project=. validation/mrgsolve/iov_cases_qspatlas.jl
Rscript validation/mrgsolve/iov_cases_mrgsolve.R && Rscript validation/mrgsolve/iov_compare.R
```

mrgsolve 2.0.1 was used. With Homebrew R on macOS, compiling mrgsolve models may need
`SDKROOT` set to an installed SDK (see the comment at the top of `cases_mrgsolve.R`).
Outputs go to `out/` next to the scripts (not committed).
