"""
    Fix(id, kind, description; evidence = "")

A documented change needed to make a model match its paper. `kind` is one of
`:typo` (equation or value misprinted), `:missing` (value or rule not stated and
inferred), `:units` (unit inconsistency), `:ambiguity` (text allows several readings;
the one used is stated), `:code` (taken from supplementary code where the text
differs), or `:from_target` (an unstated input read from, or chosen by matching, the
very figure being reproduced). `figures` names the figures the fix applies to (for
example `["fig6"]`); empty means all of them. It only affects the per-figure grades. `evidence` says where the fix comes from (supplement
file, equation number, figure) so a reader can check it; a fix without evidence
cannot earn a B.

`:from_target` fixes are circular for the values they set: a grade that depends on
one is marked with a dagger (B†) and `score.toml` records `inputs_from_target = true`.
Only use it for protocol inputs (initial conditions, which of several published
scenarios a figure shows), never for model parameters; a model whose parameters must
be read off its own figure is D.
"""
struct Fix
    id::String
    kind::Symbol
    description::String
    evidence::String
    figures::Vector{String}
    function Fix(id, kind, description; evidence = "", figures = String[])
        kind in (:typo, :missing, :units, :ambiguity, :code, :from_target) || throw(ArgumentError("unknown fix kind $kind"))
        return new(String(id), kind, String(description), String(evidence), String.(collect(figures)))
    end
end

# fixes that apply to figure `k` (a fix with no `figures` applies to all of them)
fixes_for(fixes, k) = filter(f -> isempty(f.figures) || k in f.figures, fixes)

"""
    CurveResult(ref, x, y, comparison)

A simulated curve together with its reference and comparison.
"""
struct CurveResult
    ref::ReferenceCurve
    x::Vector{Float64}
    y::Vector{Float64}
    comparison::CurveComparison
end

"""
    Run

One attempt at reproducing the figures: `label` (`"as published"` or `"with fixes"`),
`fixes` applied, `curves`, and `error` (empty if the model ran).
"""
struct Run
    label::String
    fixes::Vector{Fix}
    curves::Vector{CurveResult}
    error::String
end

ran(r::Run) = isempty(r.error)
matched(r::Run) = ran(r) && !isempty(r.curves) && all(c -> c.comparison.passed, r.curves)

"""
    run_curves(f, refs, label; fixes = Fix[], tolerance = Tolerance(), tolerances = Dict())

Build a [`Run`](@ref). `f(ref)` must return simulated `(x, y)` vectors for reference
curve `ref`; it is called once per reference. `tolerances` can override the
tolerance per curve id. Any exception inside `f` is caught and recorded as the run
error (the model did not run), so one broken model never stops the atlas build.
"""
function run_curves(f, refs, label::AbstractString; fixes = Fix[], tolerance::Tolerance = Tolerance(), tolerances = Dict{String, Tolerance}())
    curves = CurveResult[]
    try
        for ref in refs
            x, y = f(ref)
            tol = get(tolerances, ref.id, tolerance)
            push!(curves, CurveResult(ref, collect(Float64, x), collect(Float64, y), compare_curve(ref, x, y; tolerance = tol)))
        end
    catch e
        msg = sprint(showerror, e)
        return Run(String(label), collect(Fix, fixes), curves, first(msg, 2000))
    end
    return Run(String(label), collect(Fix, fixes), curves, "")
end

"""
    Reproduction

Everything known about one model's reproduction: `id`, `citation`, `doi`, `pmcid`,
the `published` run (model exactly as printed), an optional `fixed` run (after
documented [`Fix`](@ref)es), and `blocker` (why the model cannot run from published
information, for grade D; empty otherwise).
"""
Base.@kwdef struct Reproduction
    id::String
    citation::String
    doi::String
    pmcid::String = ""
    published::Union{Run, Nothing} = nothing
    fixed::Union{Run, Nothing} = nothing
    blocker::String = ""
end

"""
    grade(r::Reproduction) -> Char

- `'A'`: the model as published runs and matches every reference curve.
- `'B'`: it matches only after the documented fixes, each with evidence.
- `'C'`: it runs but does not match, with or without fixes.
- `'D'`: it cannot be run from published information (`blocker` set, or no run
  completed).

Recalibration never changes the grade: grades describe the published model, and
refitted results are reported separately.
"""
function grade(r::Reproduction)
    isempty(r.blocker) || return 'D'
    r.published !== nothing && matched(r.published) && return 'A'
    r.fixed !== nothing && matched(r.fixed) && all(f -> !isempty(strip(f.evidence)), r.fixed.fixes) && return 'B'
    any(x -> x !== nothing && ran(x), (r.published, r.fixed)) && return 'C'
    return 'D'
end

"""
    from_target(r::Reproduction) -> Bool

True when the grade rests on a `:from_target` fix (B only), i.e. some input was read
from the figure being reproduced. Shown as B† on the site.
"""
from_target(r::Reproduction) = grade(r) == 'B' && any(f -> f.kind === :from_target, r.fixed.fixes)

"""
    figure_key(ref) -> String

Which figure a reference curve belongs to: the `group` header key if present,
otherwise the figure number at the start of the curve id (`fig2_a_T` -> `fig2`,
`fig4b_neut` -> `fig4`), otherwise the id up to its first underscore.
"""
function figure_key(ref::ReferenceCurve)
    haskey(ref.meta, "group") && return ref.meta["group"]
    m = match(r"^(fig\d+)", ref.id)
    return m === nothing ? String(first(split(ref.id, '_'))) : String(m[1])
end

"""
    figure_grades(r::Reproduction) -> Vector{Pair{String, String}}

Grade per figure, using the same rules as [`grade`](@ref) on that figure's curves
alone. The model grade is set by all frozen targets together, so a model can be C
overall while one of its figures reproduces; this breakdown shows that.
"""
function figure_grades(r::Reproduction)
    isempty(r.blocker) || return Pair{String, String}[]
    runs = filter(!isnothing, [r.published, r.fixed])
    isempty(runs) && return Pair{String, String}[]
    keys_ = unique(figure_key(c.ref) for run in runs for c in run.curves)
    out = Pair{String, String}[]
    for k in keys_
        function ok(run)
            (run === nothing || !ran(run)) && return false
            cs = filter(c -> figure_key(c.ref) == k, run.curves)
            return !isempty(cs) && all(c -> c.comparison.passed, cs)   # no curves for this figure is not a pass
        end
        fk = r.fixed === nothing ? Fix[] : fixes_for(r.fixed.fixes, k)
        g = ok(r.published) ? "A" :
            (ok(r.fixed) && all(f -> !isempty(strip(f.evidence)), fk)) ?
                (any(f -> f.kind === :from_target, fk) ? "B†" : "B") : "C"
        push!(out, k => g)
    end
    return out
end

"""
    grade_label(r) -> String

`"A"`, `"B"`, `"B†"`, `"C"` or `"D"`.
"""
grade_label(r::Reproduction) = string(grade(r), from_target(r) ? "†" : "")

function run_summary(run::Union{Run, Nothing})
    run === nothing && return Dict{String, Any}()
    Dict{String, Any}(
        "label" => run.label,
        "ran" => ran(run),
        "matched" => matched(run),
        "error" => run.error,
        "fixes" => [Dict("id" => f.id, "kind" => String(f.kind), "description" => f.description, "evidence" => f.evidence,
            "figures" => f.figures) for f in run.fixes],
        "curves" => [Dict{String, Any}(
            "id" => c.ref.id, "figure" => c.ref.figure, "passed" => c.comparison.passed,
            "rms" => round(c.comparison.rms; sigdigits = 4), "p95" => round(c.comparison.p95; sigdigits = 4),
            "max" => round(c.comparison.max; sigdigits = 4),
            "tol_rms" => c.comparison.tolerance.rms, "tol_p95" => c.comparison.tolerance.p95,
            "n_points" => length(c.ref.x),
        ) for c in run.curves],
    )
end

"""
    write_score(path, r::Reproduction)

Write `score.toml` with the grade and per-curve metrics. The docs site and the test
suite read this file.
"""
function write_score(path::AbstractString, r::Reproduction)
    d = Dict{String, Any}(
        "id" => r.id, "citation" => r.citation, "doi" => r.doi, "pmcid" => r.pmcid,
        "grade" => string(grade(r)), "grade_label" => grade_label(r), "inputs_from_target" => from_target(r),
        "figure_grades" => Dict(figure_grades(r)),
        "blocker" => r.blocker,
        "published" => run_summary(r.published), "fixed" => run_summary(r.fixed),
    )
    open(io -> TOML.print(io, d; sorted = true), path, "w")
    return path
end

read_score(path::AbstractString) = TOML.parsefile(path)

function Base.show(io::IO, ::MIME"text/plain", r::Reproduction)
    println(io, r.id, ": grade ", grade_label(r), "  (", r.citation, ")")
    isempty(r.blocker) || println(io, "  blocker: ", r.blocker)
    fg = figure_grades(r)
    length(fg) > 1 && println(io, "  by figure: ", join(["$k $g" for (k, g) in fg], ", "))
    for run in (r.published, r.fixed)
        run === nothing && continue
        println(io, "  ", run.label, ": ", ran(run) ? (matched(run) ? "matches" : "runs, does not match") : "error: " * first(run.error, 200))
        for f in run.fixes
            println(io, "    fix ", f.id, " (", f.kind, "): ", f.description)
        end
        for c in run.curves
            println(io, "    ", c.comparison)
        end
    end
end
