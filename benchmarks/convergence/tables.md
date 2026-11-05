### Per model and variant: largest score change against the own settings

`d rms` and `d p95` are the largest absolute change of rms / p95 over the curves of the run (axis-box units; the pass rule is rms <= 0.015 and p95 <= 0.03 unless a curve has its own tolerance). `shift` is the largest change of any simulated point in axis-box units. `flips` counts curves whose pass/fail differs from own.

| model | run | variant | curves | max d rms | max d p95 | max shift (box) | flips | failed solves | non-finite |
|---|---|---|---|---|---|---|---|---|---|
| fornari_2019 | fixed | tol/100 | 9 | 4.9e-11 | 1.4e-10 | 2.17e-10 | 0 | 0 | 0 |
| fornari_2019 | fixed | altA Vern9 | 9 | 5e-11 | 1.4e-10 | 2.19e-10 | 0 | 0 | 0 |
| fornari_2019 | fixed | altB FBDF | 9 | 1.08e-09 | 1.92e-09 | 8.21e-09 | 0 | 0 | 0 |
| fornari_2019 | fixed | dense x4 | 9 | 4.96e-07 | 1e-06 | 0 | 0 | 0 | 0 |
| hutchinson_2016 | published | tol/100 | 30 | 2.53e-08 | 4.07e-08 | 5.04e-08 | 0 | 0 | 0 |
| hutchinson_2016 | published | altA Rodas5P | 30 | 2.78e-08 | 4.15e-08 | 5.14e-08 | 0 | 0 | 0 |
| hutchinson_2016 | published | altB Vern9 | 30 | 2.41e-08 | 3.83e-08 | 4.75e-08 | 0 | 0 | 0 |
| hutchinson_2016 | published | dense x4 | 30 | 2.32e-07 | 3.18e-07 | 0 | 0 | 0 | 0 |
| hutchinson_2016 | fixed | tol/100 | 30 | 4.25e-08 | 5.79e-08 | 6.56e-08 | 0 | 0 | 0 |
| hutchinson_2016 | fixed | altA Rodas5P | 30 | 4.25e-08 | 5.86e-08 | 7.48e-08 | 0 | 0 | 0 |
| hutchinson_2016 | fixed | altB Vern9 | 30 | 4.05e-08 | 5.52e-08 | 6.16e-08 | 0 | 0 | 0 |
| hutchinson_2016 | fixed | dense x4 | 30 | 6.68e-07 | 1.98e-06 | 0 | 0 | 0 | 0 |
| lindauer_2017 | published | tol/100 | 5 | 0 | 0 | 3.58e-09 | 0 | 3 | 0 |
| lindauer_2017 | published | altA FBDF | 5 | 0 | 0 | inf | 0 | 2 | 0 |
| lindauer_2017 | published | altB Rosenbrock23 | 5 | 1e-10 | 0 | 0.00306 | 0 | 3 | 0 |
| lindauer_2017 | published | dense x4 | 5 | 0 | 0 | 0 | 0 | 3 | 0 |
| lindauer_2017 | fixed | tol/100 | 5 | 1e-10 | 6e-10 | 1.42e-09 | 0 | 0 | 0 |
| lindauer_2017 | fixed | altA FBDF | 5 | 3.3e-09 | 1.51e-08 | 4.13e-08 | 0 | 0 | 0 |
| lindauer_2017 | fixed | altB Rosenbrock23 | 5 | 1.04e-08 | 4.71e-08 | 1.06e-07 | 0 | 0 | 0 |
| lindauer_2017 | fixed | dense x4 | 5 | 0 | 0 | 0 | 0 | 0 | 0 |
| ouerdani_2015 | published | tol/100 | 4 | 1e-11 | 1e-10 | 2.63e-10 | 0 | 0 | 0 |
| ouerdani_2015 | published | altA Rodas5P | 4 | 2e-11 | 0 | 2.42e-10 | 0 | 0 | 0 |
| ouerdani_2015 | published | altB Vern9 | 4 | 1e-11 | 1e-10 | 2.62e-10 | 0 | 0 | 0 |
| ouerdani_2015 | published | dense x4 | 4 | 1.57e-06 | 3.35e-07 | 0 | 0 | 0 | 0 |
| ouerdani_2015 | fixed | tol/100 | 4 | 4e-13 | 0 | 7.61e-11 | 0 | 0 | 0 |
| ouerdani_2015 | fixed | altA Rodas5P | 4 | 1e-11 | 1.4e-11 | 1.44e-10 | 0 | 0 | 0 |
| ouerdani_2015 | fixed | altB Vern9 | 4 | 4e-13 | 1e-12 | 7.53e-11 | 0 | 0 | 0 |
| ouerdani_2015 | fixed | dense x4 | 4 | 1.72e-07 | 1.12e-07 | 0 | 0 | 0 | 0 |
| stein_2019 | published | tol/100 | 6 | 1e-11 | 3e-11 | 2.89e-11 | 0 | 0 | 0 |
| stein_2019 | published | altA Rodas5P | 6 | 0 | 1e-11 | 2.19e-11 | 0 | 0 | 0 |
| stein_2019 | published | altB Tsit5 | 6 | 6e-11 | 1.3e-10 | 1.68e-10 | 0 | 0 | 0 |
| stein_2019 | published | dense x4 | 6 | 1e-11 | 2e-11 | 2.75e-11 | 0 | 0 | 0 |
| stein_2019 | published | reps x2 | 6 | 0.00227 | 0.00638 | 0.0103 | 0 | 0 | 0 |
| stein_2019 | published | seed+1000 | 6 | 0.00764 | 0.0105 | 0.0126 | 0 | 0 | 0 |
| stein_2019 | fixed | tol/100 | 6 | 1e-11 | 1e-11 | 5.34e-11 | 0 | 0 | 0 |
| stein_2019 | fixed | altA Rodas5P | 6 | 0 | 1e-11 | 4.14e-11 | 0 | 0 | 0 |
| stein_2019 | fixed | altB Tsit5 | 6 | 1e-11 | 1e-10 | 3.77e-10 | 0 | 0 | 0 |
| stein_2019 | fixed | dense x4 | 6 | 1e-11 | 2e-11 | 7.48e-11 | 0 | 0 | 0 |
| stein_2019 | fixed | reps x2 | 6 | 0.00345 | 0.00202 | 0.0116 | 0 | 0 | 0 |
| stein_2019 | fixed | seed+1000 | 6 | 0.00587 | 0.0042 | 0.0156 | 0 | 0 | 0 |
| stroh_2019 | published | tol/100 | 35 | 1e-10 | 5.3e-10 | 7.31e-05 | 0 | 0 | 0 |
| stroh_2019 | published | altA FBDF | 35 | 1.53e-08 | 1.59e-08 | 0.00757 | 0 | 0 | 0 |
| stroh_2019 | published | altB Rosenbrock23 | 35 | 5.06e-07 | 1.31e-06 | 0.00966 | 0 | 0 | 0 |
| stroh_2019 | published | dense x4 | 35 | 0.00254 | 0.00431 | 0 | 0 | 0 | 0 |
| stroh_2019 | fixed | tol/100 | 35 | 5.4e-11 | 5.3e-10 | 7.31e-05 | 0 | 0 | 0 |
| stroh_2019 | fixed | altA FBDF | 35 | 1.53e-08 | 1.59e-08 | 0.00757 | 0 | 0 | 0 |
| stroh_2019 | fixed | altB Rosenbrock23 | 35 | 5.06e-07 | 1.31e-06 | 0.00966 | 0 | 0 | 0 |
| stroh_2019 | fixed | dense x4 | 35 | 1.52e-07 | 2.56e-07 | 0 | 0 | 0 | 0 |

### Per curve: the three curves with the largest change of rms in each model and run

Largest change over all variants; the full per-curve, per-variant table is results.csv.

| model | run | curve | own rms | own p95 | max d rms | max d p95 | max shift (box) | variant of largest d rms |
|---|---|---|---|---|---|---|---|---|
| fornari_2019 | fixed | fig3_Plt | 0.00722 | 0.0155 | 4.96e-07 | 1e-06 | 8.21e-09 | dense x4 |
| fornari_2019 | fixed | fig4a_plt | 0.0323 | 0.054 | 3.83e-07 | 9.5e-07 | 8.53e-10 | dense x4 |
| fornari_2019 | fixed | fig4b_neut | 0.00993 | 0.0179 | 1.93e-07 | 3.1e-08 | 2.76e-09 | dense x4 |
| hutchinson_2016 | published | fig6_b_T | 0.19 | 0.263 | 2.32e-07 | 2.91e-08 | 2.93e-08 | dense x4 |
| hutchinson_2016 | published | fig2_f_V_N21 | 0.699 | 0.928 | 2.18e-07 | 3.18e-07 | 2.96e-08 | dense x4 |
| hutchinson_2016 | published | fig6_d_T | 0.231 | 0.275 | 2.06e-07 | 7.94e-08 | 2.92e-08 | dense x4 |
| hutchinson_2016 | fixed | fig2_e_T_N21 | 0.00607 | 0.0116 | 6.68e-07 | 5.94e-07 | 7.48e-08 | dense x4 |
| hutchinson_2016 | fixed | fig6_d_T | 0.244 | 0.29 | 5.87e-07 | 2.9e-09 | 3.43e-08 | dense x4 |
| hutchinson_2016 | fixed | fig2_f_V_N11 | 0.00537 | 0.00874 | 5.24e-07 | 4.09e-07 | 5.94e-08 | dense x4 |
| lindauer_2017 | published | fig4_5mgkg_x5 | 1 | 1 | 1e-10 | 0 | inf | altB Rosenbrock23 |
| lindauer_2017 | published | fig4_0p5mgkg_x0p2 | 0.939 | 1 | 0 | 0 | 4.25e-09 | tol/100 |
| lindauer_2017 | published | fig4_5mgkg_x0p2 | 0.937 | 1 | 0 | 0 | 2.03e-09 | tol/100 |
| lindauer_2017 | fixed | fig4_0p5mgkg_x5 | 0.131 | 0.359 | 1.04e-08 | 4.71e-08 | 6.28e-08 | altB Rosenbrock23 |
| lindauer_2017 | fixed | fig4_0p5mgkg_x0p2 | 0.0964 | 0.0992 | 8e-11 | 8e-11 | 4.47e-10 | altB Rosenbrock23 |
| lindauer_2017 | fixed | fig4_base_case | 0.000906 | 0.00117 | 1.74e-11 | 1e-11 | 5.51e-11 | altB Rosenbrock23 |
| ouerdani_2015 | published | fig4_right_P | 0.064 | 0.111 | 1.57e-06 | 7.5e-08 | 1.42e-10 | dense x4 |
| ouerdani_2015 | published | fig4_right_K | 0.0331 | 0.0424 | 2.4e-07 | 3.35e-07 | 4.07e-11 | dense x4 |
| ouerdani_2015 | published | fig4_left_P | 0.081 | 0.0985 | 2.11e-07 | 7.29e-08 | 2.63e-10 | dense x4 |
| ouerdani_2015 | fixed | fig4_right_K | 0.0112 | 0.016 | 1.72e-07 | 2.92e-08 | 3.96e-11 | dense x4 |
| ouerdani_2015 | fixed | fig4_right_P | 0.00882 | 0.0138 | 1.2e-07 | 8.06e-08 | 1.44e-10 | dense x4 |
| ouerdani_2015 | fixed | fig4_left_P | 0.000865 | 0.00156 | 6.71e-08 | 4.07e-08 | 7.22e-11 | dense x4 |
| stein_2019 | published | fig3a_left_p90 | 0.0441 | 0.0613 | 0.00764 | 0.0105 | 0.0119 | seed+1000 |
| stein_2019 | published | fig3a_left_p10 | 0.107 | 0.157 | 0.00298 | 0.00512 | 0.0121 | seed+1000 |
| stein_2019 | published | fig3a_right_p90 | 0.0261 | 0.0401 | 0.00192 | 0.00542 | 0.0113 | seed+1000 |
| stein_2019 | fixed | fig3a_left_p90 | 0.0322 | 0.0463 | 0.00587 | 0.0023 | 0.0129 | seed+1000 |
| stein_2019 | fixed | fig3a_left_p10 | 0.106 | 0.156 | 0.00369 | 0.00211 | 0.0061 | seed+1000 |
| stein_2019 | fixed | fig3a_left_p50 | 0.0582 | 0.0749 | 0.00204 | 0.00169 | 0.0126 | seed+1000 |
| stroh_2019 | published | fig4c_100x | 0.0573 | 0.0822 | 0.00254 | 0.00431 | 3.08e-05 | dense x4 |
| stroh_2019 | published | fig4c_10x | 0.157 | 0.178 | 0.00231 | 0.00355 | 0.00242 | dense x4 |
| stroh_2019 | published | fig4a_100x | 0.0262 | 0.0457 | 5.06e-07 | 1.31e-06 | 3.08e-05 | altB Rosenbrock23 |
| stroh_2019 | fixed | fig4a_100x | 0.0262 | 0.0457 | 5.06e-07 | 1.31e-06 | 3.08e-05 | altB Rosenbrock23 |
| stroh_2019 | fixed | fig4b_100x | 0.0307 | 0.0516 | 3.49e-07 | 1.06e-06 | 3.08e-05 | altB Rosenbrock23 |
| stroh_2019 | fixed | fig4a_10x | 0.00287 | 0.00712 | 2.01e-07 | 4.74e-08 | 0.00242 | altB Rosenbrock23 |

### Curves whose pass/fail flips in some variant, or that sit within 1.5x of a threshold

| model | run | curve | own rms | own p95 | own pass | own margin (max of rms/tol, p95/tol) | variants that flip |
|---|---|---|---|---|---|---|---|
| fornari_2019 | fixed | fig3_Neut | 0.011 | 0.0167 | true | 0.736 | none |
| ouerdani_2015 | fixed | fig4_right_K | 0.0112 | 0.016 | true | 0.745 | none |
| stroh_2019 | published | fig3a_M1S1 | 0.0192 | 0.0313 | false | 1.28 | none |
| stroh_2019 | published | fig3b_M1S2 | 0.0177 | 0.0285 | false | 1.18 | none |
| stroh_2019 | published | fig3c_M1S2 | 0.0151 | 0.0259 | false | 1.01 | none |
| stroh_2019 | published | fig3c_mAb00 | 0.00756 | 0.0235 | true | 0.784 | none |
| stroh_2019 | published | fig4b_10x | 0.0164 | 0.0212 | false | 1.1 | none |
| stroh_2019 | published | fig4d_100x | 0.0143 | 0.018 | true | 0.954 | none |
| stroh_2019 | fixed | fig3a_M1S1 | 0.0192 | 0.0313 | false | 1.28 | none |
| stroh_2019 | fixed | fig3b_M1S2 | 0.0177 | 0.0285 | false | 1.18 | none |
| stroh_2019 | fixed | fig3c_M1S2 | 0.0151 | 0.0259 | false | 1.01 | none |
| stroh_2019 | fixed | fig3c_mAb00 | 0.00756 | 0.0235 | true | 0.784 | none |
| stroh_2019 | fixed | fig4d_100x | 0.0143 | 0.018 | true | 0.954 | none |

### Solves that did not finish with a successful retcode or gave non-finite output

| model | run | solve | variant | retcode | finite |
|---|---|---|---|---|---|
| lindauer_2017 | published | 0.5:gamma:5.0 | own | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | own | Unstable | true |
| lindauer_2017 | published | 5.0:SLtg:5.0 | own | MaxIters | true |
| lindauer_2017 | published | 0.5:gamma:5.0 | tol/100 | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | tol/100 | Unstable | true |
| lindauer_2017 | published | 5.0:SLtg:5.0 | tol/100 | MaxIters | true |
| lindauer_2017 | published | 0.5:gamma:5.0 | altA FBDF | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | altA FBDF | Unstable | true |
| lindauer_2017 | published | 0.5:gamma:5.0 | altB Rosenbrock23 | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | altB Rosenbrock23 | Unstable | true |
| lindauer_2017 | published | 5.0:SLtg:5.0 | altB Rosenbrock23 | MaxIters | true |
| lindauer_2017 | published | 0.5:gamma:5.0 | dense x4 | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | dense x4 | Unstable | true |
| lindauer_2017 | published | 5.0:SLtg:5.0 | dense x4 | MaxIters | true |

### Reproduction of score.toml by the own variant

| model | curves compared | max relative rms difference | pass/fail differences |
|---|---|---|---|
| fornari_2019 | 9 | 0.000425 | 0 |
| hutchinson_2016 | 60 | 0.000372 | 0 |
| lindauer_2017 | 10 | 0.000166 | 0 |
| ouerdani_2015 | 8 | 0.000357 | 0 |
| stein_2019 | 12 | 1 | 0 |
| stroh_2019 | 70 | 0.000394 | 0 |
