using Dates
import Statistics: middle

# https://github.com/JuliaStats/Statistics.jl/issues/47
_middle(args...) = middle(args)
_middle(t1::Dates.AbstractTime) = t1
_middle(t1::T, t2::T) where {T <: Dates.AbstractTime} = T(round(Int, (t1 + t2).value / 2))

function diff_median(times)
    dts = diff(times)
    # return median!(dts)
    # copy `Statistics.median!` here to avoid type privacy
    inds = axes(dts, 1)
    n = length(inds)
    mid = div(first(inds) + last(inds), 2)
    if isodd(n)
        return _middle(partialsort!(dts, mid))
    else
        m = partialsort!(dts, mid:(mid + 1))
        return _middle(m[1], m[2])
    end
end

"""Get the time resolution of the times."""
function time_resolution(times; N = 6070)
    n = length(times)
    return n > N ? diff_median(@view(times[1:N])) : diff_median(times)
end

function time_resolution(data::DimArray)
    time_dim = dims(data, TimeDim)
    return time_resolution(parent(time_dim.val))
end

# https://docs.pola.rs/api/python/stable/reference/dataframe/api/polars.DataFrame.group_by_dynamic.html
"""
Group `x` into windows based on `every` and `period`.
"""
function groupby_dynamic(x::AbstractVector{T}, every, period = every, start_by = :window) where {T}
    min = minimum(x)
    max = maximum(x)
    group_idx = Vector{UnitRange{Int}}()
    starts = Vector{T}()
    current_start = ifelse(start_by == :window, floor(min, every), min)
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


"""
    split_range(t0, t1, n)

Split the range from `t0` to `t1` into `n` parts.
"""
function split_range(t0, t1, n::Int)
    return if n <= 1
        ((t0, t1),)
    else
        dt = (t1 - t0) / n
        ((t0 + (i - 1) * dt, min(t0 + i * dt, t1)) for i in 1:n)
    end
end

function split_range(t0, t1, dt)
    n = ceil(Int, (t1 - t0) / dt)
    return ((t0 + (i - 1) * dt, min(t0 + i * dt, t1)) for i in 1:n)
end


"""
    find_local_minimum(x, start_idx, direction, threshold)

Find local minimum starting from `start_idx` in the given `direction` (-1 for backward, +1 for forward).
First moves while above `threshold`, then continues below threshold until local minimum is found.
"""
function find_local_minimum(x, start_idx, direction, threshold)
    idx = start_idx
    min_val = x[idx]
    min_idx = idx
    n = length(x)
    
    # Helper function to check bounds
    in_bounds = direction > 0 ? <=(n) : >=(1)
    
    # Move in direction while above threshold
    while in_bounds(idx + direction) && x[idx] >= threshold
        idx += direction
        if x[idx] < min_val
            min_val = x[idx]
            min_idx = idx
        end
    end
    
    # Continue searching for local minimum below threshold
    while in_bounds(idx + direction) && x[idx] < threshold
        if x[idx] < min_val
            min_val = x[idx]
            min_idx = idx
        elseif x[idx] > min_val  # Found local minimum
            break
        end
        idx += direction
    end
    
    return min_idx
end