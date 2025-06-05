const PENALTY_FACTOR = Ref(0.1)

function _argmax_pair(metric, a)
    dist_matrix = pairwise(metric, eachrow(a))
    max_idx = argmax(dist_matrix)
    return Tuple(max_idx)
end

function argmax_pair(metric, a, ::Val{d}=Val{2}()) where {d}
    # Manually compute distances and track maximum
    max_dist = zero(eltype(a))
    max_i, max_j = 1, 1
    n = size(a, d)
    @inbounds for j in 1:n
        aj = selectdim(a, d, j)
        for i in (j+1):n
            ai = selectdim(a, d, i)
            dist = metric(ai, aj, i, j)
            if dist > max_dist
                max_dist, max_i, max_j = dist, i, j
            end
        end
    end
    return (max_i, max_j)
end

"""
    ts_max_distance(ts; query=TimeDim)

Compute the time interval when the timeseries has maximum cumulative variation.
"""
function ts_max_distance(ts; query=TimeDim, dist=Euclidean(), penalty_factor=PENALTY_FACTOR[])
    data = parent(ts)
    times = dims(ts, query)
    metric = time_penalized_metric(dist, times, penalty_factor)
    i, j = argmax_pair(metric, PermutedDimsArray(data, (2, 1)))
    return minmax(times[i], times[j])
end

"""
    time_penalized_metric(metric, times, penalty_factor)

Create a metric that penalizes long time separations with small variation.

It reduces the effective distance for long time separations.

Parameters:
- `metric`: Base distance metric function (e.g., Euclidean())
- `times`: Vector of time points
- `penalty_factor`: Controls the strength of the time penalty (higher = stronger penalty)

Returns a function that takes (ai, aj, i, j) where ai, aj are data points and i, j are indices.
"""
function time_penalized_metric(metric, times, penalty_factor)
    penalty_factor == 0 && return (ai, aj, i, j) -> metric(ai, aj)

    # Create a closure that includes time penalty in the distance calculation
    time_span = times[end] - times[1]
    return function (ai, aj, i, j)
        dist = metric(ai, aj)
        time_diff = abs(times[i] - times[j])
        normalized_time_diff = time_diff / time_span
        return dist / (1 + penalty_factor * normalized_time_diff)
    end
end