function _argmax_pair(metric, a)
    dist_matrix = pairwise(metric, eachrow(a))
    max_idx = argmax(dist_matrix)
    return Tuple(max_idx)
end

function argmax_pair(metric, a)
    # Manually compute distances and track maximum
    max_dist = -Inf
    max_i, max_j = 1, 1
    n = size(a, 2)
    @inbounds for j in 1:n
        aj = view(a, :, j)
        for i in (j+1):n
            ai = view(a, :, i)
            dist = metric(ai, aj)
            if dist > max_dist
                max_dist, max_i, max_j = dist, i, j
            end
        end
    end
    return (max_i, max_j)
end

function ts_max_distance(data, times; dist=Euclidean())
    i, j = argmax_pair(dist, PermutedDimsArray(data, (2, 1)))
    return minmax(times[i], times[j])
end

"""
    ts_max_distance(ts; query=TimeDim)

Compute the time interval when the timeseries has maximum cumulative variation.
"""
function ts_max_distance(ts; query=TimeDim, dist=Euclidean())
    data = ts.data
    times = dims(ts, query)
    return ts_max_distance(data, times; dist)
end

function process_events!(events, data; dist=Euclidean(), kwargs...)
    @chain events begin
        @rtransform! :t_us_ds = ts_max_distance(tview(data, :tstart, :tstop); dist)
    end
end