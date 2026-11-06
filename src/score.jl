"""
    Fix(id, kind, description; evidence, basis, stage, figures = String[])

A documented choice needed to run or match a model. Fields:

- `kind`: `:typo` (equation or value misprinted), `:missing` (value or rule not stated),
  `:units` (unit inconsistency), `:ambiguity` (text allows several readings; the one
  used is stated), `:code` (taken from supplementary code where the text differs), or
  `:from_target` (an unstated input read from, or chosen by matching, the very figure
  being reproduced).
- `basis`: where the value comes from. `:paper` (another part of the article or its
  tables), `:supplement` (the article's supplementary material or code), `:reference`
  (a source the article cites), `:standard` (a standard physical or physiological
  value, e.g. a molecular weight or a blood volume), `:interpretation` (a reading of
  ambiguous text), `:assumption` (a value chosen where the paper gives none or only a
  range or bound, e.g. a midpoint), `:figure` (the figure being reproduced).
  `:unspecified` is the default and is treated as outside the paper by the strict rule.
- `stage`: `:before_comparison` or `:after_comparison`, i.e. whether the choice was
  made before or after the simulation was first compared with the target. Choices
  made after comparison are reported as such.
- `evidence`: where a reader can check the value. A fix without evidence cannot earn B.
- `figures`: the figures the fix applies to (e.g. `["fig6"]`); empty means all. It
  only affects the per-figure grades.

Fixes listed in the *published* run are gap fills that run needs even before any
correction; a published run with fills can never be graded A.

`:from_target` fixes are circular for the values they set: a grade that depends on
one is marked with a dagger (B†). Only use them for protocol inputs (initial
conditions, which of several published scenarios a figure shows), never for model
parameters; a model whose parameters must be read off its own figure is D.
"""
struct Fix
    id::String
    kind::Symbol
    description::String
    evidence::String
    figures::Vector{String}
    basis::Symbol
    stage::Symbol
    function Fix(id, kind, description; evidence = "", figures = String[], basis = :unspecified, stage = :unspecified)
        kind in (:typo, :missing, :units, :ambiguity, :code, :from_target) || throw(ArgumentError("unknown fix kind $kind"))
        basis in (:paper, :supplement, :reference, :standard, :interpretation, :assumption, :figure, :unspecified) ||
            throw(ArgumentError("unknown fix basis $basis"))
        stage in (:before_comparison, :after_comparison, :unspecified) || throw(ArgumentError("unknown fix stage $stage"))
        kind === :from_target && basis ∉ (:figure, :unspecified) && throw(ArgumentError("a :from_target fix has basis :figure"))
        return new(String(id), kind, String(description), String(evidence), String.(collect(figures)), basis, stage)
    end
end

"""
    NotRunnable(msg)

Thrown by a model's simulation when the published information is not enough to run
it (an undefined symbol, a missing equation). [`run_curves`](@ref) records it as "did
not run"; any other exception is a bug in the model code and propagates.
"""
struct NotRunnable <: Exception
    msg::String
end
Base.showerror(io::IO, e::NotRunnable) = print(io, "NotRunnable: ", e.msg)

# fixes that apply to figure `k` (a fix with no `figures` applies to all of them)
fixes_for(fixes, k) = filter(f -> isempty(f.figures) || k in f.figures, fixes)

"""
    CurveResult(ref, x, y, comparison; sim_edges = nothing)

A simulated curve together with its reference, its comparison (which decides pass/fail) and
its [`CurveDiagnostics`](@ref) (report-only: coverage, reverse distance, band centre and edge
deviations). `sim_edges` is the optional `(lo = ..., hi = ...)` of [`curve_diagnostics`](@ref).
"""
struct CurveResult
    ref::ReferenceCurve
    x::Vector{Float64}
    y::Vector{Float64}
    comparison::CurveComparison
    diagnostics::CurveDiagnostics
end

CurveResult(ref, x, y, comparison; sim_edges = nothing) =
    CurveResult(ref, x, y, comparison, curve_diagnostics(ref, x, y; sim_edges))

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
tolerance per curve id. A [`NotRunnable`](@ref) thrown inside `f` is recorded as the
run error (the published information does not suffice to run the model); any other
exception propagates, so a coding error is never mistaken for a gap in the paper.

`f(ref)` may return a third element `(lo = ..., hi = ...)`, the model's own simulated interval
edges at the returned `x` (for a band reference, e.g. the `lo` and `hi` of a [`vpc`](@ref)). They
are compared with the reference edges as a report-only metric (`band_edge_rms`, `band_edge_max`
in [`CurveDiagnostics`](@ref)); the band score and the grade do not use them. Every curve's
diagnostics (coverage, reverse distance, band centre) are computed here and stored by
[`write_score`](@ref); none of them changes pass/fail. The one opt-in exception is the reverse-distance criterion of
[`Tolerance`](@ref) (`reverse_p95_max`), which is part of `compare_curve` and so of pass/fail; pass it
in `tolerance` and in every curve-specific entry of `tolerances`.
"""
function run_curves(f, refs, label::AbstractString; fixes = Fix[], tolerance::Tolerance = Tolerance(), tolerances = Dict{String, Tolerance}())
    curves = CurveResult[]
    try
        for ref in refs
            x, y, rest... = f(ref)
            sim_edges = isempty(rest) ? nothing : first(rest)
            tol = get(tolerances, ref.id, tolerance)
            push!(curves, CurveResult(ref, collect(Float64, x), collect(Float64, y), compare_curve(ref, x, y; tolerance = tol); sim_edges))
        end
    catch e
        # only a documented gap in the published information counts as "did not run";
        # anything else is a bug in the model code and must fail loudly
        e isa NotRunnable || rethrow()
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
    r.published !== nothing && matched(r.published) && isempty(r.published.fixes) && return 'A'
    matching_run(r) !== nothing && return 'B'
    any(x -> x !== nothing && ran(x), (r.published, r.fixed)) && return 'C'
    return 'D'
end

evidenced(fixes) = all(f -> !isempty(strip(f.evidence)), fixes)

# the run that earns a B: the published run if it matches with evidenced fills, else the
# fixed run if it matches with evidenced fixes
function matching_run(r::Reproduction)
    for run in (r.published, r.fixed)
        run !== nothing && matched(run) && evidenced(run.fixes) && return run
    end
    return nothing
end

const PAPER_BASES = (:paper, :supplement, :interpretation)

"""
    grade_strict(r) -> Char

The grade under the strict reading of "from published information": only the article
and its supplement count. If any run needs a value whose basis is a cited reference, a
standard value, the target figure or unspecified, the model cannot be run from the
published information alone and is D; otherwise the strict grade equals [`grade`](@ref).
"""
function grade_strict(r::Reproduction)
    g = grade(r)
    fixes = vcat([run.fixes for run in (r.published, r.fixed) if run !== nothing]...)
    any(f -> f.basis ∉ PAPER_BASES, fixes) && return 'D'
    return g
end

"""
    from_target(r::Reproduction) -> Bool

True when the grade rests on a `:from_target` fix (B only), i.e. some input was read
from the figure being reproduced. Shown as B† on the site.
"""
from_target(r::Reproduction) = grade(r) == 'B' && any(f -> f.kind === :from_target, matching_run(r).fixes)

"""
    post_comparison(r::Reproduction) -> Bool

True when the grade is B and rests on at least one fix made after the first comparison with
the target (`stage = :after_comparison`). Shown as ‡ on the site (B‡, or B†‡ together with
[`from_target`](@ref)).
"""
post_comparison(r::Reproduction) = grade(r) == 'B' && any(f -> f.stage === :after_comparison, matching_run(r).fixes)

# label for a B, from the fixes it rests on: † input read from the target, ‡ a fix made after the first comparison
b_label(fixes) = string("B", any(f -> f.kind === :from_target, fixes) ? "†" : "",
    any(f -> f.stage === :after_comparison, fixes) ? "‡" : "")

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
        # a figure is judged on its own curves: a NotRunnable raised later, for another
        # figure, does not void curves already scored
        figcurves(run) = run === nothing ? CurveResult[] : filter(c -> figure_key(c.ref) == k, run.curves)
        ok(run) = (cs = figcurves(run); !isempty(cs) && all(c -> c.comparison.passed, cs))
        if isempty(figcurves(r.published)) && isempty(figcurves(r.fixed))
            push!(out, k => "D")        # neither run produced this figure
            continue
        end
        pk = r.published === nothing ? Fix[] : fixes_for(r.published.fixes, k)
        fk = r.fixed === nothing ? Fix[] : fixes_for(r.fixed.fixes, k)
        brun = ok(r.published) && evidenced(pk) ? pk : ok(r.fixed) && evidenced(fk) ? fk : nothing
        g = ok(r.published) && isempty(pk) ? "A" :
            brun !== nothing ? b_label(brun) : "C"
        push!(out, k => g)
    end
    return out
end

"""
    grade_label(r) -> String

`"A"`, `"B"`, `"B†"`, `"B‡"`, `"B†‡"`, `"C"` or `"D"`. † marks an input read from the
figure being reproduced ([`from_target`](@ref)); ‡ marks a fix made after the first
comparison ([`post_comparison`](@ref)).
"""
grade_label(r::Reproduction) = grade(r) == 'B' ? b_label(matching_run(r).fixes) : string(grade(r))

# report-only fields of score.toml; not-applicable (NaN) values are left out
function diagnostics_summary(g::CurveDiagnostics)
    out = Dict{String, Any}("coverage" => round(g.coverage; sigdigits = 4), "sim_coverage" => round(g.sim_coverage; sigdigits = 4))
    isnan(g.fraction_checked) || (out["fraction_checked"] = round(g.fraction_checked; sigdigits = 4))
    g.rev_n > 0 && merge!(out, Dict{String, Any}("rev_rms" => round(g.rev_rms; sigdigits = 4), "rev_p95" => round(g.rev_p95; sigdigits = 4), "rev_n" => g.rev_n))
    isnan(g.band_centre_dev_rms) || merge!(out, Dict{String, Any}("band_centre_dev_rms" => round(g.band_centre_dev_rms; sigdigits = 4), "band_n" => g.band_n))
    isnan(g.band_edge_rms) || merge!(out, Dict{String, Any}("band_edge_rms" => round(g.band_edge_rms; sigdigits = 4), "band_edge_max" => round(g.band_edge_max; sigdigits = 4)))
    return out
end

# fields of the opt-in reverse-distance criterion; empty (so score.toml is unchanged) when it is off
function reverse_summary(c::CurveResult)
    tol = c.comparison.tolerance.reverse_p95_max
    tol === nothing && return Dict{String, Any}()
    out = Dict{String, Any}("tol_reverse_p95" => tol)
    c.ref.kind === :curve && (out["failed_reverse"] = c.comparison.failed_reverse)
    return out
end

function run_summary(run::Union{Run, Nothing})
    run === nothing && return Dict{String, Any}()
    Dict{String, Any}(
        "label" => run.label,
        "ran" => ran(run),
        "matched" => matched(run),
        "error" => run.error,
        "fixes" => [Dict("id" => f.id, "kind" => String(f.kind), "description" => f.description, "evidence" => f.evidence,
            "figures" => f.figures, "basis" => String(f.basis), "stage" => String(f.stage)) for f in run.fixes],
        "curves" => [Dict{String, Any}(
            "id" => c.ref.id, "figure" => c.ref.figure, "passed" => c.comparison.passed,
            "rms" => round(c.comparison.rms; sigdigits = 4), "p95" => round(c.comparison.p95; sigdigits = 4),
            "max" => round(c.comparison.max; sigdigits = 4),
            "tol_rms" => c.comparison.tolerance.rms, "tol_p95" => c.comparison.tolerance.p95,
            "n_points" => length(c.ref.x),
            diagnostics_summary(c.diagnostics)...,
            reverse_summary(c)...,
        ) for c in run.curves],
    )
end

"""
    write_score(path, r::Reproduction)

Write `score.toml` with the grade and per-curve metrics. The docs site and the test
suite read this file. Besides `rms`, `p95`, `max` and `passed` each curve carries the
report-only checks of [`CurveDiagnostics`](@ref): `coverage`, `sim_coverage`, `fraction_checked`
(not for `:bars`), `rev_rms`, `rev_p95`,
`rev_n` (`:curve` references), `band_centre_dev_rms`, `band_n`, `band_edge_rms`, `band_edge_max`
(`:band` references; the edge fields only when the model passed simulated edges). They do not
enter the grade. When a curve's [`Tolerance`](@ref) sets the opt-in `reverse_p95_max`, the entry also
has `tol_reverse_p95` and, for `:curve` references, `failed_reverse` (true when that criterion
made the curve fail, whether or not the rms and p95 rule also did).
"""
function write_score(path::AbstractString, r::Reproduction)
    d = Dict{String, Any}(
        "id" => r.id, "citation" => r.citation, "doi" => r.doi, "pmcid" => r.pmcid,
        "grade" => string(grade(r)), "grade_label" => grade_label(r), "inputs_from_target" => from_target(r),
        "grade_strict" => string(grade_strict(r)),
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
