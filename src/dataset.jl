# NONMEM-style datasets: read_nmdata, write_nmdata, NMData. Semantics and the NONMEM
# conventions followed: docs/src/simulation.md ("Datasets").

# Data items read here. Any other column is a covariate (or dropped on request).
const NM_ITEMS = (:ID, :TIME, :EVID, :AMT, :CMT, :RATE, :SS, :II, :ADDL, :MDV, :DV)
# NONMEM data items with a meaning QSPAtlas does not implement: always an error.
const NM_UNSUPPORTED = (:PCMT, :CALL, :CONT, :L1, :L2, :DATE, :DAT1, :DAT2, :DAT3,
    :XVID1, :XVID2, :XVID3, :XVID4, :XVID5, :MRG_, :RAW_, :RPT_)

"""
    NMData

A NONMEM-style dataset read by [`read_nmdata`](@ref).

  - `subjects`: one [`Subject`](@ref) per individual record, in dataset order, ready for
    [`simulate`](@ref).
  - `obs`: one row per observation record (EVID 0), in record order:
    `(id, time, endpoint, order, dv, row)`; `endpoint` is `nothing` when the record reports
    every endpoint, `order` is `:pre` or `:post`, `dv` is `missing` for ".", `row` the record number.
  - `columns`, `values`: the dataset as read (`values[j][i]` is column `j` of record `i`,
    `missing` for "." or empty), used by [`write_nmdata`](@ref).
  - `covariates`: the covariate columns; `covariate_rule`: `:nocb` or `:locf`.
"""
struct NMData
    columns::Vector{Symbol}
    values::Vector{Vector{Any}}
    subjects::Vector{Subject}
    obs::Vector{NamedTuple}
    covariates::Vector{Symbol}
    covariate_rule::Symbol
end

Base.show(io::IO, d::NMData) = print(io, "NMData(", length(d.subjects), " subjects, ", isempty(d.values) ? 0 : length(d.values[1]),
    " records, ", length(d.obs), " observations, covariates ", d.covariates, ", ", d.covariate_rule, ")")

# ---- loading ------------------------------------------------------------------------------

_missingvalue(v) = v === missing || v === nothing || (v isa AbstractString && (isempty(strip(v)) || strip(v) == "."))

function _parsefield(f::AbstractString)
    x = strip(f)
    if length(x) >= 2 && x[1] == '"' && x[end] == '"'
        x = x[2:(end - 1)]
    end
    (isempty(x) || x == ".") && return missing
    v = tryparse(Int, x)
    v === nothing || return v
    v = tryparse(Float64, x)
    v === nothing || return v
    return String(x)
end

# Columns and values (column-major) of a CSV file (comma-separated if the header has a comma,
# else whitespace-separated). The first line is the header; a leading '#' on it is dropped.
function _load_table(path::AbstractString)
    isfile(path) || throw(ArgumentError("read_nmdata: no such file $path"))
    lines = [l for l in readlines(path) if !isempty(strip(l))]
    isempty(lines) && throw(ArgumentError("read_nmdata: $path is empty"))
    head = strip(lines[1])
    startswith(head, "#") && (head = strip(head[2:end]))
    comma = occursin(',', head)
    splitline(l) = comma ? split(l, ',') : split(strip(l))
    cols = [Symbol(strip(strip(h), '"')) for h in splitline(head)]
    vals = [Any[] for _ in cols]
    for (k, l) in enumerate(lines[2:end])
        c = lstrip(l)
        (startswith(c, '#') || startswith(c, '@')) && throw(ArgumentError(
            "read_nmdata: record $k starts with '$(c[1])'; comment records (NONMEM IGNORE) are not read, remove them from the file"))
        fs = splitline(l)
        length(fs) == length(cols) || throw(ArgumentError("read_nmdata: record $k has $(length(fs)) fields, the header has $(length(cols))"))
        for j in eachindex(cols)
            push!(vals[j], _parsefield(fs[j]))
        end
    end
    return cols, vals
end

_cellvalue(v) = _missingvalue(v) ? missing : v isa AbstractString ? _parsefield(v) : v

function _load_table(rows::AbstractVector{<:NamedTuple})
    isempty(rows) && throw(ArgumentError("read_nmdata: no records"))
    cols = collect(Symbol, keys(rows[1]))
    for (k, r) in enumerate(rows)
        collect(keys(r)) == cols || throw(ArgumentError("read_nmdata: record $k has columns $(collect(keys(r))), the first record $cols"))
    end
    return cols, [Any[_cellvalue(r[c]) for r in rows] for c in cols]
end

function _load_table(tbl::NamedTuple)
    cols = collect(Symbol, keys(tbl))
    n = length(first(tbl))
    all(c -> length(c) == n, values(tbl)) || throw(ArgumentError("read_nmdata: columns have different lengths"))
    return cols, [Any[_cellvalue(v) for v in tbl[c]] for c in cols]
end

# ---- reading ------------------------------------------------------------------------------

"""
    read_nmdata(source; cmt = Dict(), endpoints = nothing, covariates = nothing, drop = Symbol[],
                covariate_rule = :nocb, control = DoseControl(), occasion = nothing) -> NMData

Read a NONMEM-style dataset: a CSV file path (comma- or whitespace-separated, header first,
"." or empty for missing), a vector of NamedTuples (one per record) or a NamedTuple of columns.
Returns an [`NMData`](@ref) whose `subjects` keep the record order of the file.

Data items (upper or lower case): `ID`, `TIME` (required), `EVID` (0 observation, 1 dose,
2 other, 3 reset, 4 reset and dose), `AMT`, `CMT`, `RATE` (0 bolus, > 0 rate, -1 modelled rate,
-2 modelled duration), `SS` (0, 1, 2), `II`, `ADDL`, `MDV`, `DV`. Without `EVID` it is derived as
the NONMEM help describes: a record with `AMT`, `RATE` or `SS` nonzero is a dose (so is a
steady-state infusion, `AMT` 0, `RATE` > 0, `SS` 1); any other record is an other-type record
if `MDV` is 1, else an observation. Every other column is a covariate
(numeric, never missing), unless listed in `drop`; with `covariates` given, any other column
is an error. Unsupported items (`PCMT`, `CALL`, `CONT`, `L1`, `L2`, `DATE`, ...) and values
(EVID > 4, SS = 3, a decreasing TIME, an observation with MDV = 1, an ID that reappears after
another ID, ...) are errors, never ignored.

  - `cmt`: dose compartment of each `CMT` value, `Dict(1 => :depot, 2 => :central)`. A text
    `CMT` without an entry is used as the route name itself; a number without one is an error.
  - `endpoints`: endpoint of each observation `CMT` value, `Dict(2 => :conc)`, or `:all` to report
    every endpoint at every observation. Observations with `CMT` 0 or missing report every endpoint.
  - `covariate_rule`: `:nocb` (default, NONMEM: a record's covariate values apply over the
    interval ending at that record) or `:locf` (they apply from the record onward). Records at
    the same time must agree on every covariate.
  - `control`: a [`DoseControl`](@ref), or a function `(id, covariates) -> DoseControl` (the
    baseline covariates), for bioavailability, lag and the modelled rates and durations that
    `RATE = -1 / -2` need. Dose control that depends on random effects, on time-varying
    covariates or on the occasion belongs in the model (`PopulationModel(...; dosecontrol)`),
    which evaluates it on each dose record with that record's covariates.
  - `occasion`: the name of the occasion column (for example `:OCC`), for between-occasion
    variability. Its values are integer occasion labels, never missing; records at one time
    must agree. The occasion of a record follows `covariate_rule` like a covariate (under NOCB
    the occasion on a record applies over the interval ending at it), and a dose uses the
    occasion on its record. Without `occasion` an OCC column is an ordinary covariate.

Record order is kept: an observation record placed before a dose, reset or steady-state record
at the same time sees the state before it (`:pre`), one placed after it the state after it.
Doses implied by `ADDL` at the time of an observation record come before it. An observation
between two such records at one time, a reset or SS=1 record after a dose record at the same
time, and an observation before a dose record at a time where an `ADDL` dose also falls are errors.
"""
function read_nmdata(source; cmt = Dict(), endpoints = nothing, covariates = nothing, drop = Symbol[],
        covariate_rule::Symbol = :nocb, control = DoseControl(), occasion = nothing)
    covariate_rule in (:nocb, :locf) || throw(ArgumentError("covariate_rule must be :nocb or :locf, got $(repr(covariate_rule))"))
    endpoints === nothing || endpoints === :all || endpoints isa AbstractDict ||
        throw(ArgumentError("endpoints must be nothing, :all or a Dict, got $(repr(endpoints))"))
    cols, vals = _load_table(source)
    allunique(cols) || throw(ArgumentError("read_nmdata: duplicated column names $cols"))
    n = isempty(vals) ? 0 : length(vals[1])
    n == 0 && throw(ArgumentError("read_nmdata: no records"))

    # classify columns
    item = Dict{Symbol, Int}()
    covcols = Symbol[]
    dropped = Set(Symbol.(drop))
    occcol = nothing
    if occasion !== nothing
        occcol = findfirst(c -> c == Symbol(occasion), cols)
        occcol === nothing && throw(ArgumentError("read_nmdata: occasion column $(occasion) is not a column"))
        Symbol(uppercase(String(cols[occcol]))) in NM_ITEMS && throw(ArgumentError("read_nmdata: occasion column $(occasion) is a NONMEM data item"))
    end
    for (j, c) in enumerate(cols)
        u = Symbol(uppercase(String(c)))
        if j === occcol
            continue
        elseif u in NM_ITEMS
            haskey(item, u) && throw(ArgumentError("read_nmdata: data item $u appears twice (columns $(cols[item[u]]) and $c)"))
            item[u] = j
        elseif u in NM_UNSUPPORTED
            throw(ArgumentError("read_nmdata: data item $c is not supported"))
        elseif c in dropped
            continue
        elseif covariates === nothing || c in Symbol.(covariates)
            push!(covcols, c)
        else
            throw(ArgumentError("read_nmdata: column $c is not a NONMEM data item read here and not in `covariates`; list it in `covariates` or `drop`"))
        end
    end
    if covariates !== nothing
        for c in Symbol.(covariates)
            c in covcols || throw(ArgumentError("read_nmdata: covariate $c is not a column"))
        end
    end
    for req in (:ID, :TIME)
        haskey(item, req) || throw(ArgumentError("read_nmdata: the dataset has no $req column"))
    end
    raw(u, i) = haskey(item, u) ? vals[item[u]][i] : missing

    recs = [_nmrecord(i, raw, item, vals, cols, covcols, cmt, endpoints) for i in 1:n]
    occs = occcol === nothing ? nothing : map(1:n) do i
        v = vals[occcol][i]
        (v isa Real && isfinite(v) && isinteger(v)) || throw(ArgumentError(
            "read_nmdata: record $i: occasion $(cols[occcol]) must be an integer, got $(v === missing ? "missing" : repr(v))"))
        Int(v)
    end

    # individual records: contiguous IDs
    groups = Vector{Vector{Int}}()
    seen = Set{String}()
    for i in 1:n
        if i == 1 || recs[i].id != recs[i - 1].id
            recs[i].id in seen && throw(ArgumentError(
                "read_nmdata: ID $(recs[i].id) reappears at record $i after other IDs (NONMEM would make it a new individual); give it a different ID"))
            push!(seen, recs[i].id)
            push!(groups, Int[])
        end
        push!(groups[end], i)
    end
    subjects = Subject[]
    obs = NamedTuple[]
    for g in groups
        s, o = _nmsubject(recs[g], covcols, covariate_rule, control, occs === nothing ? nothing : occs[g])
        push!(subjects, s)
        append!(obs, o)
    end
    return NMData(cols, vals, subjects, obs, covcols, covariate_rule)
end

function _num(v, c, i)
    v === missing && return 0.0                                   # NONMEM reads "." and blank as 0
    v isa Real && isfinite(v) && return Float64(v)
    throw(ArgumentError("read_nmdata: record $i, $c: expected a finite number, got $(repr(v))"))
end
function _int(v, c, i)
    x = _num(v, c, i)
    isinteger(x) || throw(ArgumentError("read_nmdata: record $i, $c: expected an integer, got $(repr(v))"))
    return Int(x)
end

# CMT value -> Symbol via a user mapping: `nothing` if 0 or missing, an error if a number has no entry.
function _cmtname(v, mapping, what, i)
    v === missing && return nothing
    if v isa Real
        x = _int(v, "CMT", i)
        x == 0 && return nothing
        x < 0 && throw(ArgumentError("read_nmdata: record $i: negative CMT ($x, turn a compartment off) is not supported"))
        for (k, s) in (mapping === nothing ? () : mapping)
            k isa Real && k == x && return Symbol(s)
        end
        throw(ArgumentError("read_nmdata: record $i: CMT $x has no $what; pass $(what == "route" ? "cmt" : "endpoints") = Dict($x => :name)"))
    end
    name = String(v)
    for (k, s) in (mapping === nothing ? () : mapping)
        (k isa AbstractString || k isa Symbol) && String(k) == name && return Symbol(s)
    end
    return Symbol(name)
end

_idstring(v) = v isa Real && isinteger(v) ? string(Int(v)) : string(v)

# One record, checked and normalised.
function _nmrecord(i, raw, item, vals, cols, covcols, cmtmap, endpoints)
    idv = raw(:ID, i)
    idv === missing && throw(ArgumentError("read_nmdata: record $i has no ID"))
    tv = raw(:TIME, i)
    tv === missing && throw(ArgumentError("read_nmdata: record $i has no TIME"))
    time = _num(tv, "TIME", i)
    amt = _num(raw(:AMT, i), "AMT", i)
    rate = _num(raw(:RATE, i), "RATE", i)
    ss = _int(raw(:SS, i), "SS", i)
    ii = _num(raw(:II, i), "II", i)
    addl = _int(raw(:ADDL, i), "ADDL", i)
    mdvv = raw(:MDV, i)
    mdv = haskey(item, :MDV) ? _int(mdvv, "MDV", i) : -1
    mdv in (-1, 0, 1) || throw(ArgumentError("read_nmdata: record $i: MDV must be 0 or 1, got $mdv"))
    # Without EVID (NONMEM help, EVID): a record with AMT, RATE or SS nonzero is a dose event, so a
    # steady-state infusion (AMT 0, RATE > 0, SS 1) is one too; any other record is an
    # observation, or an other-type record (EVID 2) when MDV is 1
    evid = if haskey(item, :EVID)
        raw(:EVID, i) === missing && throw(ArgumentError("read_nmdata: record $i has no EVID"))
        _int(raw(:EVID, i), "EVID", i)
    elseif amt != 0 || rate != 0 || ss != 0
        1
    else
        mdv == 1 ? 2 : 0
    end
    evid in 0:4 || throw(ArgumentError("read_nmdata: record $i: EVID $evid is not supported (0 to 4)"))
    # MDV: 0 on observation records, 1 on every other record (NONMEM)
    if mdv >= 0
        evid == 0 && mdv == 1 && throw(ArgumentError(
            "read_nmdata: record $i: an observation record (EVID 0) with MDV 1 is not supported (NONMEM requires MDV 0 on observation records); use EVID 2 for a record without an observation"))
        evid != 0 && mdv == 0 && throw(ArgumentError(
            "read_nmdata: record $i: EVID $evid with MDV 0 is not supported (MDV must be 1 on non-observation records)"))
    end
    ss in (0, 1, 2) || throw(ArgumentError("read_nmdata: record $i: SS $ss is not supported (0, 1 or 2)"))
    ii >= 0 || throw(ArgumentError("read_nmdata: record $i: II must be >= 0, got $ii"))
    addl >= 0 || throw(ArgumentError("read_nmdata: record $i: ADDL must be >= 0, got $addl"))
    isdose = evid == 1 || evid == 4
    if !isdose
        for (nm, v) in (("AMT", amt), ("RATE", rate), ("SS", ss), ("II", ii), ("ADDL", addl))
            v == 0 || throw(ArgumentError("read_nmdata: record $i: $nm = $v on an EVID $evid record (must be 0 or missing)"))
        end
    end
    cmtv = raw(:CMT, i)
    route = nothing
    endpoint = nothing
    if isdose
        route = _cmtname(cmtv, cmtmap, "route", i)
        route === nothing && throw(ArgumentError("read_nmdata: record $i: a dose record needs a CMT (non-zero)"))
        rate in (-1.0, -2.0) || rate >= 0 || throw(ArgumentError("read_nmdata: record $i: RATE $rate is not supported (0, > 0, -1 or -2)"))
        evid == 4 && ss != 0 && throw(ArgumentError("read_nmdata: record $i: SS = $ss on a reset-and-dose record (EVID 4) is not supported"))
        constant = ss != 0 && amt == 0 && (rate > 0 || rate == -1)
        if constant
            ii == 0 || throw(ArgumentError("read_nmdata: record $i: a steady-state infusion (AMT 0, RATE $rate) needs II = 0"))
            addl == 0 || throw(ArgumentError("read_nmdata: record $i: ADDL must be 0 on a steady-state infusion record"))
        else
            amt > 0 || throw(ArgumentError("read_nmdata: record $i: a dose record needs AMT > 0 (AMT 0 is only allowed for a steady-state infusion: SS 1 or 2, RATE > 0 or -1, II 0)"))
            ss != 0 && ii <= 0 && throw(ArgumentError("read_nmdata: record $i: a steady-state dose (SS $ss) needs II > 0"))
            addl > 0 && ii <= 0 && throw(ArgumentError("read_nmdata: record $i: ADDL $addl needs II > 0"))
            ii > 0 && ss == 0 && addl == 0 && throw(ArgumentError("read_nmdata: record $i: II $ii on a dose without SS or ADDL has no meaning; set it to 0"))
        end
    elseif evid == 0
        endpoint = endpoints === :all ? nothing : _cmtname(cmtv, endpoints, "endpoint", i)
    elseif evid == 2
        cmtv isa Real && cmtv < 0 && throw(ArgumentError("read_nmdata: record $i: EVID 2 with a negative CMT (turn a compartment off) is not supported"))
        cmtv isa AbstractString && throw(ArgumentError("read_nmdata: record $i: EVID 2 with CMT $(repr(cmtv)) is not supported (compartments are always on; use 0 or leave it empty)"))
    end
    dvv = raw(:DV, i)
    dv = evid == 0 && dvv !== missing ? _num(dvv, "DV", i) : missing
    cov = NamedTuple{Tuple(covcols)}(Tuple(begin
        v = vals[findfirst(==(c), cols)][i]
        (v === missing || !(v isa Real)) && throw(ArgumentError(
            "read_nmdata: record $i: covariate $c is $(v === missing ? "missing" : repr(v)); covariates must be numbers (NONMEM would read a missing value as 0)"))
        Float64(v)
    end for c in covcols))
    return (row = i, id = _idstring(idv), time, evid, amt, rate, ss, ii, addl, route, endpoint, dv, cov)
end

# One individual record -> Subject and its observation rows.
function _nmsubject(recs, covcols, rule, control, occs = nothing)
    id = recs[1].id
    for k in 2:length(recs)
        recs[k].time >= recs[k - 1].time || throw(ArgumentError(
            "read_nmdata: ID $id, record $(recs[k].row): TIME $(recs[k].time) is before the previous record's $(recs[k - 1].time) (time may not decrease, also not after a reset)"))
    end
    isevent(r) = r.evid in (1, 3, 4)

    # same-time order that the regimen cannot express
    for r in recs
        (r.evid in (3, 4) || r.ss == 1) || continue
        for q in recs
            q.row < r.row && q.time == r.time && q.evid in (1, 4) && throw(ArgumentError(
                "read_nmdata: ID $id, record $(r.row): a $(r.evid == 1 ? "SS=1 dose" : "reset") record after a dose record at the same time (record $(q.row)) is not supported (NONMEM would cancel that dose; QSPAtlas applies resets and steady states first)"))
        end
    end

    # regimen, in record order
    pieces = Any[]
    for r in recs
        if r.evid == 1 || r.evid == 4
            rate = r.rate == -1 ? :model : r.rate > 0 ? r.rate : 0.0
            duration = r.rate == -2 ? :model : 0.0
            push!(pieces, Regimen(r.amt; time = r.time, cmt = r.route, ii = r.ii, addl = r.addl, rate, duration, ss = r.ss, reset = r.evid == 4))
        elseif r.evid == 3
            push!(pieces, Reset(r.time))
        end
    end
    regimen = try
        isempty(pieces) ? Regimen(Dose[]) : Regimen(pieces...)
    catch e
        e isa ArgumentError ? throw(ArgumentError("read_nmdata: ID $id: $(e.msg)")) : rethrow()
    end
    implied = Set{Float64}()                      # nominal times of ADDL doses
    for r in recs, k in 1:r.addl
        push!(implied, r.time + k * r.ii)          # same arithmetic as Regimen
    end

    # observations, tagged by record order
    ot = Float64[]
    oo = Symbol[]
    oe = Union{Nothing, Symbol}[]
    rows = NamedTuple[]
    for r in recs
        r.evid == 0 || continue
        before = any(q -> isevent(q) && q.time == r.time && q.row < r.row, recs)
        after = any(q -> isevent(q) && q.time == r.time && q.row > r.row, recs)
        before && after && throw(ArgumentError(
            "read_nmdata: ID $id, record $(r.row): an observation between two dose or reset records at the same time is not supported"))
        after && r.time in implied && throw(ArgumentError(
            "read_nmdata: ID $id, record $(r.row): an observation before a dose record at time $(r.time), where an ADDL dose also falls, is not supported (the ADDL dose comes before the observation, the record after it)"))
        order = after ? :pre : :post
        push!(ot, r.time)
        push!(oo, order)
        push!(oe, r.endpoint)
        push!(rows, (id = id, time = r.time, endpoint = r.endpoint, order, dv = r.dv, row = r.row))
    end

    # covariates: one value per distinct time
    times = Float64[]
    cvals = NamedTuple[]
    for r in recs
        if !isempty(times) && times[end] == r.time
            r.cov == cvals[end] || throw(ArgumentError(
                "read_nmdata: ID $id, record $(r.row): covariates differ from another record at the same time $(r.time) (not supported; records at one time must agree)"))
        else
            push!(times, r.time)
            push!(cvals, r.cov)
        end
    end
    tvcols = [c for c in covcols if length(unique(v[c] for v in cvals)) > 1]
    tv = nothing
    if !isempty(tvcols)
        sub(v) = NamedTuple{Tuple(tvcols)}(Tuple(v[c] for c in tvcols))
        m = length(times)
        keep = rule === :nocb ? [k == m || sub(cvals[k]) != sub(cvals[k + 1]) for k in 1:m] :
               [k > 1 && sub(cvals[k]) != sub(cvals[k - 1]) for k in 1:m]
        tv = [times[k] => sub(cvals[k]) for k in 1:m if keep[k]]
    end
    # occasions: one label per distinct time, kept where it changes (as for covariates)
    occasions = nothing
    if occs !== nothing
        otimes = Float64[]
        olab = Int[]
        for (r, o) in zip(recs, occs)
            if !isempty(otimes) && otimes[end] == r.time
                o == olab[end] || throw(ArgumentError(
                    "read_nmdata: ID $id, record $(r.row): occasion $o differs from another record at the same time $(r.time) ($(olab[end]))"))
            else
                push!(otimes, r.time)
                push!(olab, o)
            end
        end
        m = length(otimes)
        keep = rule === :nocb ? [k == m || olab[k] != olab[k + 1] for k in 1:m] : [k == 1 || olab[k] != olab[k - 1] for k in 1:m]
        occasions = [otimes[k] => olab[k] for k in 1:m if keep[k]]
    end
    base = cvals[1]
    ctl = control isa DoseControl ? control : control(id, base)
    ctl isa DoseControl || throw(ArgumentError("read_nmdata: control must return a DoseControl, got $(typeof(ctl))"))
    s = Subject(id; regimen, covariates = base, obstimes = ot, control = ctl, timevarying = tv,
        obsorder = oo, endpoints = oe, covariate_rule = rule, occasions)
    return s, rows
end

# ---- writing ------------------------------------------------------------------------------

_csvfield(v) = v === missing ? "." : v isa AbstractString ? (occursin(r"[,\"\s]", v) ? "\"" * replace(v, "\"" => "\"\"") * "\"" : v) :
    v isa AbstractFloat && isinteger(v) && abs(v) < 1e15 ? string(Int(v)) : string(v)

"""
    write_nmdata(path, data::NMData; result = nothing, replicate = 1, value = :dv, endpoint = nothing)

Write the dataset as a comma-separated file with the original columns, "." for missing values.
With a [`SimResult`](@ref) from simulating `data.subjects`, the `DV` column of every observation
record is replaced by the simulated `value` (`:dv` or `:ipred`) of replicate `replicate`: record
by record, matching each subject's observation records in order. Observation records that
report every endpoint need `endpoint` to say which one goes into `DV`.
"""
function write_nmdata(path::AbstractString, data::NMData; result = nothing, replicate::Integer = 1, value::Symbol = :dv, endpoint = nothing)
    vals = [copy(v) for v in data.values]
    if result !== nothing
        result isa SimResult || throw(ArgumentError("result must be a SimResult"))
        value in (:dv, :ipred) || throw(ArgumentError("value must be :dv or :ipred"))
        j = findfirst(c -> uppercase(String(c)) == "DV", data.columns)
        j === nothing && throw(ArgumentError("write_nmdata: the dataset has no DV column"))
        sims = Dict{Tuple{String, Symbol}, Vector{Float64}}()
        for o in result.obs
            o.replicate == replicate && push!(get!(sims, (o.id, o.endpoint), Float64[]), Float64(getproperty(o, value)))
        end
        used = Dict{Tuple{String, Symbol}, Int}()
        for o in data.obs
            ep = o.endpoint === nothing ? endpoint : o.endpoint
            ep === nothing && throw(ArgumentError("write_nmdata: record $(o.row) reports every endpoint; pass `endpoint` to choose the one written to DV"))
            key = (o.id, Symbol(ep))
            k = get(used, key, 0) + 1
            v = get(sims, key, Float64[])
            k <= length(v) || throw(ArgumentError("write_nmdata: the result has no value for record $(o.row) (ID $(o.id), endpoint $ep, replicate $replicate)"))
            used[key] = k
            vals[j][o.row] = v[k]
        end
    end
    open(path, "w") do io
        println(io, join(String.(data.columns), ","))
        for i in eachindex(vals[1])
            println(io, join((_csvfield(vals[j][i]) for j in eachindex(vals)), ","))
        end
    end
    return path
end
