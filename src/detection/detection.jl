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

"""
    detect_variance(f, tmin, tmax, tau; split=1, kwargs...)

Apply variance detection on a time series `f(tmin, tmax)`

Optional splitting into multiple chunks to improve performance and memory efficiency.
"""
function detect_variance(f, tmin, tmax, tau; split=nothing, kwargs...)
    if isnothing(split)
        return detect_variance(f(tmin, tmax; add_unit=false), tau; kwargs...)
    else
        tranges = split_range(tmin, tmax, split)
        @showprogress mapreduce(append!, tranges) do (chunk_min, chunk_max)
            detect_variance(f, chunk_min, chunk_max, tau; kwargs...)
        end
    end
end

function ids_finder(f, tmin, tmax, tau; split=nothing, kwargs...)
    if isnothing(split)
        data = f(tmin, tmax; add_unit=false)
        events = detect_variance(data, tau; kwargs...)
        process_events!(events, data; kwargs...)
    else
        tranges = split_range(tmin, tmax, split)
        @showprogress mapreduce(append!, tranges) do (chunk_min, chunk_max)
            ids_finder(f, chunk_min, chunk_max, tau; kwargs...)
        end
    end
end