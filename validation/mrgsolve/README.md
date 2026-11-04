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

Run from the repository root:

```bash
julia --startup-file=no --threads=4 --project=. validation/mrgsolve/cases_qspatlas.jl
julia --startup-file=no --threads=4 --project=. validation/mrgsolve/pop_qspatlas.jl
Rscript validation/mrgsolve/cases_mrgsolve.R && Rscript validation/mrgsolve/compare.R
Rscript validation/mrgsolve/pop_mrgsolve.R && Rscript validation/mrgsolve/pop_compare.R
```

mrgsolve 2.0.1 was used. With Homebrew R on macOS, compiling mrgsolve models may need
`SDKROOT` set to an installed SDK (see the comment at the top of `cases_mrgsolve.R`).
Outputs go to `out/` next to the scripts (not committed).
