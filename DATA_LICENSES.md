# Data licences and attribution

Status: draft record prepared 2026-11-05. It is a factual inventory, not legal advice, and
it makes no legal conclusion. Open points are listed at the end.

## Scope

- The MIT licence in `LICENSE` applies to the code in this repository (the Julia package in
  `src/`, tests, tools, scripts, and the original code in `models/*/`) and to the original
  text and plots that we wrote.
- It does not relicense any published article, figure, table or supplement.
- Digitised coordinates (`models/*/reference/*.csv`, `digitise.toml`), transcribed tables and
  matrices, and parameter values restated from the papers are facts extracted from the cited
  articles. They are not our own creation, and they remain subject to the articles' terms
  where those terms apply. The MIT licence does not widen what the article licences allow.
- Paper figure images and supplement files are not in the repository. The digitising scripts
  read them from a local folder (`tools/digitise/README.md`).
- Each article's licence below was read from the `<license>` element of its PubMed Central
  XML file (`papers/<id>/PMC*.xml`). The XML carries one licence statement per article. It
  does not state a separate licence for the supplementary files.

## Summary

| model | article licence | licence URL | copyright line in the XML |
|---|---|---|---|
| ouerdani_2015 | CC BY-NC-ND 4.0 | https://creativecommons.org/licenses/by-nc-nd/4.0/ | (c) 2015 The Authors, CPT: Pharmacometrics & Systems Pharmacology published by Wiley Periodicals, Inc. on behalf of the American Society for Clinical Pharmacology and Therapeutics |
| lindauer_2017 | CC BY-NC-ND 4.0 | http://creativecommons.org/licenses/by-nc-nd/4.0/ | (c) 2016 The Authors, same publisher statement |
| stroh_2019 | CC BY-NC-ND 4.0 | https://creativecommons.org/licenses/by-nc-nd/4.0/ | (c) 2019 CytomX Therapeutics, Inc., same publisher statement |
| stein_2019 | CC BY-NC-ND 4.0 | https://creativecommons.org/licenses/by-nc-nd/4.0/ | (c) 2019 The Authors, same publisher statement |
| fornari_2019 | CC BY-NC 4.0 | http://creativecommons.org/licenses/by-nc/4.0/ | (c) 2019 Astrazeneca (spelling as in the XML), same publisher statement |
| hutchinson_2016 | CC BY 4.0 | https://creativecommons.org/licenses/by/4.0/ | (c) 2016 The Authors, same publisher statement |

The licence text in each XML (quoted):

- ouerdani_2015: "This is an open access article under the terms of the Creative Commons
  Attribution-NonCommercial-NoDerivs License, which permits use and distribution in any
  medium, provided the original work is properly cited, the use is non-commercial and no
  modifications or adaptations are made." (`license-type="creativeCommonsBy-nc-nd"`)
- lindauer_2017: same sentence as ouerdani_2015 (no `license-type` attribute; link
  http://creativecommons.org/licenses/by-nc-nd/4.0/).
- stroh_2019: "This is an open access article under the terms of the
  http://creativecommons.org/licenses/by-nc-nd/4.0/ License, which permits use and
  distribution in any medium, provided the original work is properly cited, the use is
  non-commercial and no modifications or adaptations are made."
  (`license-type="creativeCommonsBy-nc-nd"`). The XML also lists "John Wiley & Sons, Ltd" as
  `<copyright-holder>` on other elements (not the article-level statement).
- stein_2019: same sentence as stroh_2019 (`license-type="creativeCommonsBy-nc-nd"`; Wiley is
  also named as `<copyright-holder>` on other elements).
- fornari_2019: "This is an open access article under the terms of the
  http://creativecommons.org/licenses/by-nc/4.0/ License, which permits use, distribution and
  reproduction in any medium, provided the original work is properly cited and is not used for
  commercial purposes."
- hutchinson_2016: "This is an open access article under the terms of the Creative Commons
  Attribution License, which permits use, distribution and reproduction in any medium,
  provided the original work is properly cited." (`license-type="creativeCommonsBy"`)

The four NC-ND articles and the NC article carry the non-commercial condition. Anyone who
reuses the files listed below should check that condition against their own use. Using a
public MIT repository does not by itself remove it.

## Per-paper records

### ouerdani_2015

| item | value |
|---|---|
| Article | Ouerdani A, Struemper H, Suttle AB, Ouellet D, Ribba B. Preclinical Modeling of Tumor Growth and Angiogenesis Inhibition to Describe Pazopanib Clinical Effects in Renal Cell Carcinoma. CPT Pharmacometrics Syst Pharmacol 2015;4(11):660-668. DOI 10.1002/psp4.12001. PMC4716582 |
| Licence | CC BY-NC-ND 4.0, https://creativecommons.org/licenses/by-nc-nd/4.0/ |
| Copyright | (c) 2015 The Authors, published by Wiley Periodicals, Inc. on behalf of ASCPT |
| Used from the paper | Equations and parameter values (main text, tables, supplement NONMEM code); Figure 4 (typical tumour curves, mouse and patient) as digitised coordinates |
| Repo files derived from it | `models/ouerdani_2015/reference/fig4_left_K.csv`, `fig4_left_P.csv`, `fig4_right_K.csv`, `fig4_right_P.csv`; `reference/digitise.toml` (pixel coordinates of every picked point and calibration); `model.jl`, `notes.md`, `report.md` (parameter values and equations restated); `figures/fig4.png` (our regenerated Figure 4); `figures/comparison.png` (our simulation against the digitised points) |
| Also touched | `examples/population_ouerdani.jl` and `examples/population_ouerdani.png` use the paper's parameter values and variability; `benchmarks/solvers/results/ouerdani_2015*` are our solver outputs on this model |

### lindauer_2017

| item | value |
|---|---|
| Article | Lindauer A, Valiathan CR, Mehta K, Sriram V, de Greef R, Elassaiss-Schaap J, de Alwis DP. Translational Pharmacokinetic/Pharmacodynamic Modeling of Tumor Growth Inhibition Supports Dose-Range Selection of the Anti-PD-1 Antibody Pembrolizumab. CPT Pharmacometrics Syst Pharmacol 2017;6(1):11-20. DOI 10.1002/psp4.12130. PMC5270293 |
| Licence | CC BY-NC-ND 4.0, http://creativecommons.org/licenses/by-nc-nd/4.0/ |
| Copyright | (c) 2016 The Authors, published by Wiley Periodicals, Inc. on behalf of ASCPT |
| Used from the paper | Equations and parameters (Table 1, supplement control streams); Figure 4 (tornado plot, 64 bars) as digitised bar tips; two values read from the Figure 4 caption text |
| Repo files derived from it | `models/lindauer_2017/reference/fig4_0p5mgkg_x0p2.csv`, `fig4_0p5mgkg_x5.csv`, `fig4_5mgkg_x0p2.csv`, `fig4_5mgkg_x5.csv` (bar tips); `fig4_base_case.csv` (values from the caption text, not digitised); `reference/digitise.toml` (hand-picked pixel positions); `model.jl`, `notes.md`, `report.md` (parameters and equations restated); `figures/figure4_regenerated.png` (our redraw of the tornado plot); `figures/comparison.png` |
| Not copied | The mouse dataset in the supplement (`s007.csv`) is not in the repository (stated in `notes.md`) |

### stroh_2019

| item | value |
|---|---|
| Article | Stroh M, Sagert J, Burke JM, Apgar JF, Lin L, Millard BL, Kavanaugh WM. Quantitative Systems Pharmacology Model of a Masked, Tumor-Activated Antibody. CPT Pharmacometrics Syst Pharmacol 2019;8(9):676-684. DOI 10.1002/psp4.12448. PMC6765697 |
| Licence | CC BY-NC-ND 4.0, https://creativecommons.org/licenses/by-nc-nd/4.0/ |
| Copyright | (c) 2019 CytomX Therapeutics, Inc., published by Wiley Periodicals, Inc. on behalf of ASCPT |
| Used from the paper | Equations and parameters (Tables 1 and 2); Figures 3a-c and 4a-e as digitised coordinates; the supplement model code (`PSP4-8-676-s003.pdf`, KroneckerBio matrices A1 and A2 and a 32-entry parameter list), transcribed |
| Repo files derived from it | 35 CSVs in `models/stroh_2019/reference/` (`fig3a_*`, `fig3b_*`, `fig3c_*`, `fig4a_*` to `fig4e_*`); `reference/digitise.toml` (pixel boxes, tick calibrations, colours; generated by `make_toml.jl`); `model.jl`, `notes.md`, `report.md` (parameters, reaction rules read from the supplement matrices, restated); `figures/fig3_regenerated.png`, `figures/fig4_regenerated.png`, `figures/comparison.png` |
| Transcribed supplement data | `models/stroh_2019/supplement_check/A1.csv` (234 nonzero entries of A1), `A2.csv` (84 nonzero entries of A2), `k.csv` (the 32 printed parameter values); `check_matrices.jl` compares them with `model.jl` |

### stein_2019

| item | value |
|---|---|
| Article | Stein AM, Grupp SA, Levine JE, Laetsch TW, Pulsipher MA, Boyer MW, August KJ, Levine BL, Tomassian L, Shah S, Leung M, Huang P-H, Awasthi R, Mueller KT, Wood PA, June CH. Tisagenlecleucel Model-Based Cellular Kinetic Analysis of Chimeric Antigen Receptor-T Cells. CPT Pharmacometrics Syst Pharmacol 2019;8(5):285-295. DOI 10.1002/psp4.12388. PMC6539725 |
| Licence | CC BY-NC-ND 4.0, https://creativecommons.org/licenses/by-nc-nd/4.0/ |
| Copyright | (c) 2019 The Authors, published by Wiley Periodicals, Inc. on behalf of ASCPT |
| Used from the paper | Equations and parameters (text, tables, supplement Monolix code); Figure 3a (shaded simulated percentile bands, left and right panels) as digitised band edges |
| Repo files derived from it | `models/stein_2019/reference/fig3a_left_p10.csv`, `fig3a_left_p50.csv`, `fig3a_left_p90.csv`, `fig3a_right_p10.csv`, `fig3a_right_p50.csv`, `fig3a_right_p90.csv`; `reference/digitise.toml`; `model.jl`, `notes.md`, `report.md` (parameters and equations restated); `figures/fig3.png` (our simulated bands); `figures/comparison.png`, `figures/comparison_fixed.png` |
| Not copied | The supplementary data set (`PSP4-8-285-s002.csv`) is read from a local folder by `reproduce.jl` if present and is not in the repository |

### fornari_2019

| item | value |
|---|---|
| Article | Fornari C, Oplustil O'Connor L, Pin C, Smith A, Yates JWT, Cheung SYA, Jodrell DI, Mettetal JT, Collins TA. Quantifying Drug-Induced Bone Marrow Toxicity Using a Novel Haematopoiesis Systems Pharmacology Model. CPT Pharmacometrics Syst Pharmacol 2019;8(11):858-868. DOI 10.1002/psp4.12459. PMC6875710 |
| Licence | CC BY-NC 4.0, http://creativecommons.org/licenses/by-nc/4.0/ |
| Copyright | (c) 2019 Astrazeneca (as written in the XML), published by Wiley Periodicals, Inc. on behalf of ASCPT |
| Used from the paper | Equations and parameters (Table 1, supplement tables S1 to S4 and model code); Figures 3 and 4 as digitised coordinates; Supplement Figure S3 (rat free platinum and Damage curves) as digitised coordinates |
| Repo files derived from it | `models/fornari_2019/reference/fig3_CMP.csv`, `fig3_MEP.csv`, `fig3_MPP.csv`, `fig3_Neut.csv`, `fig3_Plt.csv`, `fig4a_neut.csv`, `fig4a_plt.csv`, `fig4b_neut.csv`, `fig4b_plt.csv`; `reference/ancillary/S3_damage.csv`, `S3_freePt_uM.csv` (from the supplement figure); `reference/digitise.toml`; `dose_scale_from_S3.jl` (diagnostic that fits a scale to the S3 curves and restates Table S3 values); `pk_units_check.jl`; `model.jl`, `notes.md`, `report.md`, `sensitivity.jl` (parameters restated); `figures/fig3_regenerated.png`, `figures/fig4_regenerated.png`, `figures/comparison.png` |

### hutchinson_2016

| item | value |
|---|---|
| Article | Hutchinson LG, Mueller H-J, Gaffney EA, Maini PK, Wagg J, Phipps A, Boetsch C, Byrne HM, Ribba B. Modeling Longitudinal Preclinical Tumor Size Data to Identify Transient Dynamics in Tumor Response to Antiangiogenic Drugs. CPT Pharmacometrics Syst Pharmacol 2016;5(11):636-645. DOI 10.1002/psp4.12142. PMC5192995 |
| Licence | CC BY 4.0, https://creativecommons.org/licenses/by/4.0/ |
| Copyright | (c) 2016 The Authors, published by Wiley Periodicals, Inc. on behalf of ASCPT |
| Used from the paper | Equations and parameters (text, supplement MLXTRAN code); Figure 2 (panels a to f) and Figure 6 (panels a to d) as digitised curves; Figure 1d (group mean tumour volume markers) as digitised points for recalibration |
| Repo files derived from it | 30 CSVs in `models/hutchinson_2016/reference/` (`fig2_a_T.csv`, `fig2_b_V.csv`, `fig2_c_*`, `fig2_d_*`, `fig2_e_*`, `fig2_f_*`, `fig6_a_*` to `fig6_d_*`); `reference/digitise.toml`, `reference/digitise_fig6.toml`; `benchmarks/calibration/data/hutchinson_fig1d.csv` and `benchmarks/calibration/digitise_fig1d.jl`; `model.jl`, `notes.md`, `report.md` (parameters and equations restated); `figures/fig2.png`, `figures/fig6.png` (our regenerated figures); `figures/comparison.png`; outputs in `benchmarks/calibration/results/` (fits to the digitised Figure 1d means) |

## Other tracked files that mention the papers

No paper content is in these beyond citations and restated parameter values:
`benchmarks/solvers/results/*` (our solver timings and plots), `benchmarks/calibration/results/*`
(our fits, see hutchinson_2016), `docs/`, `tools/digitise/` (the digitiser; its README states
that images are not distributed), `validation/`, `test/`.
`models/README.md` and the top level `README.md` show example reference-CSV headers and
results tables.

## Attribution statements

The following statements are meant to be kept with any redistribution of the listed files.
Changes made by us, common to all: values were read from the published figures as pixel
coordinates and converted to axis units by calibration (about 1 pixel of reading error, see
each CSV header); the models were re-implemented in Julia from the published equations and
tables; plots were regenerated from our simulations.

**Ouerdani 2015.** Digitised coordinates and restated parameters derived from: A. Ouerdani,
H. Struemper, A. B. Suttle, D. Ouellet, B. Ribba, "Preclinical Modeling of Tumor Growth and
Angiogenesis Inhibition to Describe Pazopanib Clinical Effects in Renal Cell Carcinoma",
CPT: Pharmacometrics & Systems Pharmacology 4(11), 660-668 (2015),
https://doi.org/10.1002/psp4.12001. (c) 2015 The Authors, licensed under CC BY-NC-ND 4.0,
https://creativecommons.org/licenses/by-nc-nd/4.0/. Changes: coordinates of Figure 4 were
digitised from the figure and are not the authors' numerical data; parameters were restated.

**Lindauer 2017.** Digitised bar tips and restated parameters derived from: A. Lindauer,
C. R. Valiathan, K. Mehta, V. Sriram, R. de Greef, J. Elassaiss-Schaap, D. P. de Alwis,
"Translational Pharmacokinetic/Pharmacodynamic Modeling of Tumor Growth Inhibition Supports
Dose-Range Selection of the Anti-PD-1 Antibody Pembrolizumab", CPT: Pharmacometrics &
Systems Pharmacology 6(1), 11-20 (2017), https://doi.org/10.1002/psp4.12130. (c) 2016 The
Authors, licensed under CC BY-NC-ND 4.0, http://creativecommons.org/licenses/by-nc-nd/4.0/.
Changes: bar tips of Figure 4 were digitised by hand from the figure; two values were read
from its caption; parameters were restated.

**Stroh 2019.** Digitised coordinates, restated parameters and transcribed supplement
matrices derived from: M. Stroh, J. Sagert, J. M. Burke, J. F. Apgar, L. Lin, B. L. Millard,
W. M. Kavanaugh, "Quantitative Systems Pharmacology Model of a Masked, Tumor-Activated
Antibody", CPT: Pharmacometrics & Systems Pharmacology 8(9), 676-684 (2019),
https://doi.org/10.1002/psp4.12448. (c) 2019 CytomX Therapeutics, Inc., licensed under
CC BY-NC-ND 4.0, https://creativecommons.org/licenses/by-nc-nd/4.0/. Changes: Figures 3 and 4
were digitised; parameters were restated; the nonzero entries of the supplement matrices A1 and
A2 and the printed parameter list were transcribed from `PSP4-8-676-s003.pdf` into CSV
(see Open questions).

**Stein 2019.** Digitised percentile bands and restated parameters derived from: A. M. Stein,
S. A. Grupp, J. E. Levine, T. W. Laetsch, M. A. Pulsipher, M. W. Boyer, K. J. August,
B. L. Levine, L. Tomassian, S. Shah, M. Leung, P.-H. Huang, R. Awasthi, K. T. Mueller,
P. A. Wood, C. H. June, "Tisagenlecleucel Model-Based Cellular Kinetic Analysis of Chimeric
Antigen Receptor-T Cells", CPT: Pharmacometrics & Systems Pharmacology 8(5), 285-295 (2019),
https://doi.org/10.1002/psp4.12388. (c) 2019 The Authors, licensed under CC BY-NC-ND 4.0,
https://creativecommons.org/licenses/by-nc-nd/4.0/. Changes: band edges of Figure 3a were
digitised; parameters were restated.

**Fornari 2019.** Digitised coordinates and restated parameters derived from: C. Fornari,
L. Oplustil O'Connor, C. Pin, A. Smith, J. W. T. Yates, S. Y. A. Cheung, D. I. Jodrell,
J. T. Mettetal, T. A. Collins, "Quantifying Drug-Induced Bone Marrow Toxicity Using a Novel
Haematopoiesis Systems Pharmacology Model", CPT: Pharmacometrics & Systems Pharmacology 8(11),
858-868 (2019), https://doi.org/10.1002/psp4.12459. (c) 2019 Astrazeneca, licensed under
CC BY-NC 4.0, http://creativecommons.org/licenses/by-nc/4.0/. Changes: Figures 3 and 4 and
Supplement Figure S3 were digitised; parameters from the article and its supplement were
restated, and some were re-expressed (for example per-kg PK values, see the model report).

**Hutchinson 2016.** Digitised coordinates and restated parameters derived from: L. G.
Hutchinson, H.-J. Mueller, E. A. Gaffney, P. K. Maini, J. Wagg, A. Phipps, C. Boetsch,
H. M. Byrne, B. Ribba, "Modeling Longitudinal Preclinical Tumor Size Data to Identify
Transient Dynamics in Tumor Response to Antiangiogenic Drugs", CPT: Pharmacometrics & Systems
Pharmacology 5(11), 636-645 (2016), https://doi.org/10.1002/psp4.12142. (c) 2016 The Authors,
licensed under CC BY 4.0, https://creativecommons.org/licenses/by/4.0/. Changes: Figures 1d,
2 and 6 were digitised; parameters were restated.

## Open questions

These are unresolved. None is answered here.

1. **Digitised coordinates of NC-ND figures.** For ouerdani_2015, lindauer_2017, stroh_2019
   and stein_2019 the licence text says "no modifications or adaptations are made". It is
   not settled whether coordinates read off a figure (curves, band edges, bar tips) are
   adapted material of the figure, or numerical facts about the published result. The
   reference CSVs are what the grading uses, so this affects the core of the repository for
   four of six models.
2. **Restated parameters and equations.** Parameter tables and equations are restated in
   `model.jl`, `notes.md` and `report.md`. The same question applies to them, although they
   are shorter and cited line by line.
3. **Regenerated figures.** `figures/*_regenerated.png`, `fig2.png`, `fig3.png`, `fig4.png`,
   `fig6.png` redraw the paper's figures from our simulations in the paper's layout. It is
   open whether a redraw in the same layout is an adaptation for the NC-ND papers. The
   overlay images that sit on the original figures (`figures/digitise_check*.png`) are
   excluded by `.gitignore` for this reason.
4. **Non-commercial condition and the MIT licence.** The MIT licence permits commercial use.
   The CSVs for five of six papers are derived from NC or NC-ND material. A downstream user
   reading only `LICENSE` could miss this. Options: keep this file prominent (done in the
   README), or move the reference data out of the MIT scope explicitly.
5. **Supplementary files.** The XML states one licence per article and none for the
   supplements. Fornari 2019 Supplement Figure S3 (`reference/ancillary/S3_*.csv`) and the
   Stroh 2019 supplement matrices come from supplements. Whether the article licence covers
   them is not shown by the XML.
6. **Stroh 2019 supplement matrices.** `supplement_check/A1.csv`, `A2.csv` and `k.csv` are
   entry-by-entry transcriptions of the supplement `PSP4-8-676-s003.pdf`. Two conflicts:
   - They sit awkwardly with the ND term (a transcription of a supplement is arguably a
     copy or adaptation of it), more clearly than the figure coordinates do.
   - The README said "no ... supplementary files are redistributed here". The matrices are
     not the PDF, but they carry its numerical content, so the sentence was not accurate
     as written (the README paragraph has since been corrected).

   Recommendation: do not keep them in the public tree. The preferred option is to replace
   them with a script (for example `supplement_check/extract_matrices.jl`) that reads
   `$QSPATLAS_PAPERS/stroh_2019/supplement/PSP4-8-676-s003.pdf` from the user's own
   download and writes the three CSVs into an ignored local folder, with `check_matrices.jl`
   reading from there and skipping with a message when the files are absent. Record the
   result of the check (234 of 234 A1 entries and 84 of 84 A2 entries identical) in
   `notes.md` and `report.md`, which is already done. Add the three paths to `.gitignore`.
   Second choice: remove them and keep only the check result. Keeping them is the weakest
   option, defensible only if the right holders or the publisher confirm that matrices from
   the supplement may be redistributed. The files are in the first public commit, so
   removing them from the working tree does not remove them from history. If that matters,
   a history rewrite would be needed. The files have not been deleted yet; this is a
   decision for the repository owner.
7. **Hutchinson 2016 (CC BY).** No restriction beyond attribution applies. The attribution
   statement above should be kept with the CSVs and with `benchmarks/calibration/data/`.
8. **Fornari 2019 (CC BY-NC).** Adaptations are allowed with attribution and indication of
   changes. Only the non-commercial condition applies (see question 4).
9. **Copyright line for Fornari 2019.** The XML reads "Astrazeneca". The attribution
   statement keeps that spelling. It should be confirmed against the journal page.
10. **Per-file attribution.** CC licences ask for attribution with the material. The CSV headers carry
    the source figure and article identifier. They do not each carry the licence line. A
    `reference/README.md` per model now points here. Whether that is enough is open.
11. **Contact.** The README says the readings have not been checked with the authors. No
    permission has been requested from authors or publisher for any item above. Asking the
    corresponding authors of the four NC-ND papers is a possible way to settle questions 1
    to 3.
