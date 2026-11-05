# Figure digitiser

Every reference curve in this repository (`models/<id>/reference/*.csv`) was read off a figure of the
published article with the tools in this folder. The settings that produced each CSV (pixel box,
tick calibration, colour, excluded regions) are stored next to it in `digitise.toml`, so any reference
can be regenerated and checked.

## What it does

`digitise.jl` loads a figure image, keeps the pixels of one colour inside a pixel box for one panel,
takes the median row per pixel column, converts pixels to data values with a two-tick calibration per
axis (linear or log), thins the result and writes a CSV in the format read by
`QSPAtlas.load_reference` (header lines `# id`, `# figure`, `# x`, `# y`, `# kind`, `# source`,
`# method`, then `x,y`). It also has helpers for finding calibration pixels in a REPL
(`axes_lines`, `colour_census`, `column_profile`, `row_profile`).

`digitise.jl` can be run on any `digitise.toml` that has `[[curve]]` tables
(`julia --project=tools/digitise tools/digitise/digitise.jl <toml>`). Each model also has a script in
`models/<id>/reference/` that adds what its figures need (rectangle exclusions for legends, hue-based
separation of overlapping lines, manual rules for black solid and dashed lines, bar-tip picking); they
use the same calibration conventions and write the same CSV format. `papers.jl` holds the shared
lookup of the image folder.

## The figure images are not in this repository

The articles are published under Creative Commons licences (some CC BY-NC-ND, which does not allow
adaptations), so to stay clear of licence questions the repository contains only the digitised numbers
and the scripts, not the figures. To regenerate a reference you need the images locally.

Scripts read images from

    $QSPATLAS_PAPERS/<model id>/supplement/<file>

where `QSPATLAS_PAPERS` is an environment variable. If it is not set, the default is a folder
`papers` next to the repository folder (`joinpath(repo_root, "..", "papers")`). A script that cannot
find its image stops with a message naming the missing path.

## Getting the images

Europe PMC serves the supplementary files of an open-access article as a zip that includes the figure
JPGs:

    https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/supplementaryFiles

For example

    mkdir -p papers/ouerdani_2015/supplement
    curl -L -o ouerdani.zip https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4716582/supplementaryFiles
    unzip -j ouerdani.zip -d papers/ouerdani_2015/supplement

| model | PMCID | image files used |
|---|---|---|
| `ouerdani_2015` | PMC4716582 | `PSP4-4-660-g004.jpg` |
| `hutchinson_2016` | PMC5192995 | `PSP4-5-636-g002.jpg`, `PSP4-5-636-g006.jpg` (and `PSP4-5-636-g001.jpg` for `benchmarks/calibration`) |
| `lindauer_2017` | PMC5270293 | `PSP4-6-11-g004.jpg` (bar tips are stored in `digitise.toml`, the image is only needed to re-check them) |
| `fornari_2019` | PMC6875710 | `PSP4-8-858-g003.jpg`, `PSP4-8-858-g004.jpg`, and `PSP4-8-858-s003.pdf` (Figure S3, needs `pdftoppm`) |
| `stein_2019` | PMC6539725 | `PSP4-8-285-g003.jpg` (optional local diagnostic: `PSP4-8-285-s002.csv`, used by `models/stein_2019/reproduce.jl` only if present) |
| `stroh_2019` | PMC6765697 | `PSP4-8-676-g003.jpg`, `PSP4-8-676-g004.jpg` |

Check the licence of each article before redistributing anything you download.

## Image provenance keys

Every table in a `digitise*.toml` that names an image with `image = ...` also carries two keys next to
it, so the exact file behind a reference can be found and checked:

| key | meaning |
|---|---|
| `image_url` | public Europe PMC URL the image was obtained from: the `supplementaryFiles` zip of the article (the URL shown above). The image is the zip member named by `image` (for `fornari_2019` the short names `g003`/`g004` mean `PSP4-8-858-g003.jpg`/`-g004.jpg`). |
| `image_sha256` | SHA-256 of the exact file used for the digitisation, as hex |

The hash was checked against a fresh download of each zip (2026-11-05): the members are byte-identical
to the local copies the references were traced from. To verify a local copy:

    shasum -a 256 $QSPATLAS_PAPERS/<model id>/supplement/<file>     # compare with image_sha256

The tools ignore these keys (they read only the keys they need), and the generated files keep them:
`models/ouerdani_2015/reference/trace_fig4.jl` writes them into its `digitise.toml` and
`models/stroh_2019/reference/make_toml.jl` into each curve, so regenerating gives the same files.
The Figure 1d image used by `benchmarks/calibration/digitise_fig1d.jl` has no toml; its URL and hash are in the
header of that script. The Fornari Figure S3 PDF (`PSP4-8-858-s003.pdf`, ancillary curves only) has
SHA-256 `f67179e46047f659653255702a1f94560e8814e31e09ec41e7e9158a4de05fac` from the same zip.

## Regenerating a model's references

Run from the repository root with the `tools/digitise` environment (first time:
`julia --project=tools/digitise -e 'using Pkg; Pkg.instantiate()'`). The scripts write into the
model's `reference/` folder, so `git diff` shows whether the regenerated CSVs differ from the
committed ones.

    # one script per model (all of them include tools/digitise/digitise.jl and tools/digitise/papers.jl)
    julia --project=tools/digitise models/ouerdani_2015/reference/trace_fig4.jl [output folder]
    julia --project=tools/digitise models/stein_2019/reference/trace_fig3a.jl
    julia --project=tools/digitise models/hutchinson_2016/reference/digitise_hutchinson.jl models/hutchinson_2016/reference/digitise.toml
    julia --project=tools/digitise models/hutchinson_2016/reference/digitise_hutchinson.jl models/hutchinson_2016/reference/digitise_fig6.toml
    julia --project=tools/digitise models/stroh_2019/reference/digitise_stroh.jl
    julia --project=tools/digitise models/fornari_2019/reference/digitise_run.jl
    julia --project=tools/digitise models/lindauer_2017/reference/digitise_bars.jl models/lindauer_2017/reference/digitise.toml
    julia --project=tools/digitise benchmarks/calibration/digitise_fig1d.jl

`models/<id>/reference/check_overlay.jl` (and `digitise_check.jl` / `check_digitisation.jl` for
`stroh_2019` and `hutchinson_2016`) draw the digitised points over the original figure and read back
calibration ticks that were not used for the calibration. The overlay images contain the papers'
figures, so they are written to `models/<id>/figures/digitise_check.png` and are not committed.
