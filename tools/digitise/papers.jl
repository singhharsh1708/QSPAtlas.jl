# Where the figure images live. The images are NOT part of this repository (the articles are
# licensed CC BY-NC or CC BY-NC-ND), so every digitising script reads them from a local folder:
#
#   <QSPATLAS_PAPERS>/<model id>/supplement/<file>
#
# QSPATLAS_PAPERS is an environment variable; if unset it defaults to a `papers` folder next to
# the repository folder (i.e. joinpath(repo_root, "..", "papers")). See tools/digitise/README.md.

"""Root folder with the downloaded figure images (`ENV["QSPATLAS_PAPERS"]`, default `<repo>/../papers`)."""
papers_root() = abspath(expanduser(get(ENV, "QSPATLAS_PAPERS", joinpath(@__DIR__, "..", "..", "..", "papers"))))

"""Folder with the images of one model: `<papers_root>/<id>/supplement`."""
paper_dir(id) = joinpath(papers_root(), id, "supplement")

"""Full path of one image file of model `id`; stops with a pointer to the README if it is missing."""
function paper_image(id, file)
    path = joinpath(paper_dir(id), file)
    isfile(path) || error("figure image not found: $path\n" *
        "Figure images are not distributed with this repository. Download them as described in " *
        "tools/digitise/README.md and set QSPATLAS_PAPERS (currently $(papers_root())).")
    return path
end

"""Image folder named by a digitise.toml: `paper = \"<id>\"` (preferred), else `image_root`, else the toml's own folder."""
image_root(cfg, tomlpath) = haskey(cfg, "paper") ? paper_dir(cfg["paper"]) :
    expanduser(get(cfg, "image_root", dirname(tomlpath)))
