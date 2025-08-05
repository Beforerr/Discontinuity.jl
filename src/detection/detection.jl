using Statistics, NaNStatistics
using StaticArrays
using OhMyThreads
using ProgressMeter: @showprogress
using Distances: pairwise, Euclidean
using SPEDAS.TimeseriesUtilities: other_dims, times

export detect_variance, PVI
export process_events!
export ids_finder

include("utils.jl")
include("PVI.jl")
include("variance.jl")
include("features.jl")

"""
    ids_finder(data, tau, [V, n]; kwargs...)

Find the discontinuities in the time series `data` with the given `tau`.

Provide additional arguments `V` and `n` to calculate the features of the events.
"""
function ids_finder(data::AbstractArray, tau, args...; process = NamedTuple(), kwargs...)
    events = detect_variance(data, tau; kwargs...)
    determine_boundary!(events, data)
    return process_events!(events, data, args...; process...)
end

function ids_finder_PVI(data::AbstractArray, args...; process = NamedTuple(), kwargs...)
    PVI_ts = PVI(data; kwargs...)
    cs_intervals = coherent_intervals(PVI_ts; gap = 1)
    df = DataFrame(cs_intervals)
    determine_boundary!(df, data)
    select_current_sheets!(df, data; stability_threshold = 0.3, consistency_threshold = 0.8)
    return process_events!(df, data, args...; process...)
end

for f in (:detect_variance,)

    doc = """
        $(f)(f, tmin, tmax, tau; split=1, kwargs...)

    Apply `$(f)` on a time series `f(tmin, tmax)`

    Optional splitting into multiple chunks to improve performance and memory efficiency.
    """

    @eval @doc $doc function $f(f, tmin, tmax, tau; split = nothing, kwargs...)
        return if isnothing(split)
            $f(f(tmin, tmax; add_unit = false), tau; clean = true, kwargs...)
        else
            tranges = split_range(tmin, tmax, split)
            @showprogress mapreduce(append!, tranges) do trange
                $f(f, trange..., tau; clean = true, kwargs...)
            end
        end
    end
end


"""
    ids_finder(f, tmin, tmax, tau, args...; split=1, kwargs...)

Apply `ids_finder` on a time series `f(tmin, tmax)`

Optional splitting into multiple chunks to improve performance and memory efficiency.
"""
function ids_finder(fB, tmin, tmax, tau, args...; split = nothing, kwargs...)
    return if isnothing(split)
        new_args = map(f -> f(tmin, tmax), args)
        data = fB(tmin, tmax; add_unit = false)
        isnothing(data) && return DataFrame()
        ids_finder(data, tau, new_args...; clean = true, kwargs...)
    else
        tranges = split_range(tmin, tmax, split)
        @showprogress mapreduce(append!, tranges) do trange
            ids_finder(fB, trange..., tau, args...; kwargs...)
        end
    end
end
