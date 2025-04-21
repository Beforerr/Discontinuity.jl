using SPEDAS: mva, mva_eigen, rotate
using PlasmaFormulary: inertial_length

include("./features/mva.jl")

const SV3 = SVector{3}

function _argmax_pair(metric, a)
    dist_matrix = pairwise(metric, eachrow(a))
    max_idx = argmax(dist_matrix)
    return Tuple(max_idx)
end

function argmax_pair(metric, a::AbstractMatrix{T}) where T
    # Manually compute distances and track maximum
    max_dist = zero(T)
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
    data = parent(ts)
    times = dims(ts, query)
    return ts_max_distance(data, times; dist)
end

function stat_features(data)
    B_mags = map(norm, eachrow(data))
    B_mag = nanmean(B_mags)
    dBmag_over_Bmag = abs(B_mags[end] - B_mags[1]) / B_mag
    return (; B_mag, dBmag_over_Bmag)
end

function process_events!(events, data; dist=Euclidean(), B_unit=DEFAULT_B_UNIT, kwargs...)
    @chain events begin
        @rtransform! :t_us_ds = ts_max_distance(tview(data, :tstart, :tstop); dist)
        @rtransform! @astable begin
            t_us = :t_us_ds[1]
            t_ds = :t_us_ds[2]
            data_us = SV3(tview(data, t_us))
            data_ds = SV3(tview(data, t_ds))
            :duration = t_ds - t_us
            :n_cross = cross(data_us, data_ds)
            :ω = angle_between(data_us, data_ds)
        end
        @rtransform! $AsTable = stat_features(tview(data, :t_us_ds...))
        @rtransform! $AsTable = mva_features(tview(data, :t_us_ds...))
        transform!(
            [:B_mag, :B_lmn_before, :B_lmn_after] .=> unitize(B_unit),
            renamecols=false
        )
        @transform! :θ_mva_cross = angle_between.(:n_mva, :n_cross)
    end
end

function process_events!(events, data, V; kwargs...)
    @chain begin
        process_events!(events, data; kwargs...)
        @rtransform! @astable begin
            V_t = tview(V, :time)
            :V = V_t
            # t_us = :t_us_ds[1]
            # t_ds = :t_us_ds[2]
            # :V_us = tview(tview(V, :tstart, t_us), t_us)
            # :V_ds = tview(tview(V, t_ds, :tstop), t_ds)
        end
        @transform! @astable begin
            :V_n_cross = sproj.(:V, :n_cross)
            :V_n_mva = sproj.(:V, :n_mva)
            :L_n_cross = :duration .* :V_n_cross
            :L_n_mva = :duration .* :V_n_mva
            :J_m_max_mva = gradient_current.(:grad, :V_n_mva)
            :J_m_max_cross = gradient_current.(:grad, :V_n_cross)
        end
    end
end


function process_events!(events, data, V, n; kwargs...)
    @chain begin
        process_events!(events, data, V; kwargs...)
        @rtransform! @astable begin
            n_t = only(tview(n, :time))
            :n = n_t
        end
        @transform! @astable begin
            :d_i = inertial_length.(:n, Unitful.q, Unitful.mp)
            :V_A = Alfven_speed.(:B_mag, :n) # Alfven speed
            :J_A = @. upreferred(:V_A * :n * Unitful.q)
            :V_A_lmn_before = alfven_velocity.(:B_lmn_before, :n)
            :V_A_lmn_after = alfven_velocity.(:B_lmn_after, :n)
        end
    end
end
