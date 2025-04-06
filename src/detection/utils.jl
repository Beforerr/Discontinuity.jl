"""
    time_resolution(data; dim=Ti)

Get the time resolution of the data.
"""
function time_resolution(data; dim=Ti)
    time_dim = dims(data, dim)
    times = parent(time_dim.val)
    dts = diff(times)
    dt0 = oneunit(eltype(dts))
    return median(dts ./ dt0) * dt0
end

function groupby_dynamic(x::AbstractVector{T}, every, period=every) where T
    min = minimum(x)
    max = maximum(x)
    group_idx = Vector{UnitRange{Int}}()
    starts = Vector{T}()
    current_start = min
    while current_start <= max
        window_end = current_start + period
        # Find indices of rows that fall in the current window using searchsorted for better performance
        start_idx = searchsortedfirst(x, current_start)
        end_idx = searchsortedlast(x, window_end)
        if start_idx <= end_idx
            indices = start_idx:end_idx
            push!(group_idx, indices)
            push!(starts, current_start)
        end
        current_start += every
    end
    return group_idx, starts
end

groupby_dynamic(x::Dimension, args...; kwargs...) =
    groupby_dynamic(parent(x.val), args...; kwargs...)
