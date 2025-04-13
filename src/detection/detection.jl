using DataFrames, DataFramesMeta
using DimensionalData
using DimensionalData: dims, TimeDim
using DimensionalData.Dimensions: Dimension
using LinearAlgebra
using Statistics
using NaNStatistics
using StaticArrays
using OhMyThreads
using ProgressMeter
using Distances: pairwise, Euclidean

export detect_variance
export process_events!
export ids_finder

include("utils.jl")
include("variance.jl")
include("features.jl")

function ids_finder(data, tau; kwargs...)
    events = detect_variance(data, tau; kwargs...)
    process_events!(events, data; kwargs...)
end

for f in (:detect_variance, :ids_finder)

    doc = """
        $(f)(f, tmin, tmax, tau; split=1, kwargs...)

    Apply `$(f)` on a time series `f(tmin, tmax)`

    Optional splitting into multiple chunks to improve performance and memory efficiency.
    """

    @eval @doc $doc function $f(f, tmin, tmax, tau; split=nothing, kwargs...)
        if isnothing(split)
            $f(f(tmin, tmax; add_unit=false), tau; kwargs...)
        else
            tranges = split_range(tmin, tmax, split)
            @showprogress mapreduce(append!, tranges) do trange
                $f(f, trange..., tau; kwargs...)
            end
        end
    end
end