function _argmax_pair(metric, a)
    dist_matrix = pairwise(metric, eachrow(a))
    max_idx = argmax(dist_matrix)
    return Tuple(max_idx)
end

function argmax_pair(metric, a, ::Val{d}=Val{2}()) where {
    d}
    # Manually compute distances and track maximum
    max_dist = zero(eltype(a))
    max_i, max_j = 1, 1
    n = size(a, d)
    @inbounds for j in 1:n
        aj = selectdim(a, d, j)
        for i in (j+1):n
            ai = selectdim(a, d, i)
            dist = metric(ai, aj)
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
function ts_max_distance(ts; query=TimeDim, dist=Euclidean())
    data = parent(ts)
    times = dims(ts, query)
    i, j = argmax_pair(dist, PermutedDimsArray(data, (2, 1)))
    return minmax(times[i], times[j])
end