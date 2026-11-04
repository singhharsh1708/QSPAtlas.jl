using DiffEqCallbacks: PresetTimeCallback
using SymbolicIndexingInterface: getu, setu

"""
    repeated_doses(start, interval, n) -> Vector{Float64}

Dose times `start, start + interval, ...` (`n` doses).
"""
repeated_doses(start::Real, interval::Real, n::Integer) = [Float64(start + i * interval) for i in 0:(n - 1)]

"""
    bolus_callback(sys, state, times, amounts)

Callback that adds `amounts[i]` (or a scalar amount) to `state` at each of `times`.
`state` is a symbolic unknown of the compiled system `sys`, for example
`sys.depot`. Use it for IV bolus or for oral doses into a depot compartment.

A dose at the initial time is applied by the callback too, so leave the initial
condition of `state` at its pre-dose value. Doses at the same time all apply (they add up).
"""
function bolus_callback(sys, state, times, amounts)
    ts = collect(Float64, times)
    amts = amounts isa Real ? fill(Float64(amounts), length(ts)) : collect(Float64, amounts)
    length(amts) == length(ts) || throw(DimensionMismatch("times and amounts differ in length"))
    # Doses at the same time all apply (see dosing_callback).
    events = [DoseEvent(ts[i], :bolus, :bolus, amts[i], 0.0) for i in eachindex(ts)]
    cb, _ = dosing_callback(sys, events, [Route(:bolus, state)])
    return cb
end
