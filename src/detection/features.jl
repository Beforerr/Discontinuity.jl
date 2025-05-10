using SPEDAS: mva, mva_eigen, rotate

const SV3 = SVector{3}

include("./features/duration.jl")
include("./features/mva.jl")

function stat_features(data)
    B_mags = map(norm, eachrow(data))
    B_mag = nanmean(B_mags)
    dBmag_over_Bmag = abs(B_mags[end] - B_mags[1]) / B_mag
    return (; B_mag, dBmag_over_Bmag)
end

@views function cross_features(data)
    data_us = SV3(data[1, :])
    data_ds = SV3(data[end, :])
    n_cross = normalize(cross(data_us, data_ds))
    B_n_cross = sproj(nanmean(data, dims=1), n_cross)
    ω = angle_between(data_us, data_ds)
    return (; n_cross, B_n_cross, ω)
end

function process_events!(events, data; dist=Euclidean(), kwargs...)
    @chain events begin
        @rtransform! :t_us_ds = ts_max_distance(tview(data, :tstart, :tstop); dist)
        @rtransform! $AsTable = cross_features(parent(tview(data, :t_us_ds...)))
        @rtransform! $AsTable = stat_features(tview(data, :t_us_ds...))
        @rtransform! $AsTable = mva_features(tview(data, :t_us_ds...))
    end
end

# transform!(
# [:B_mag, :B_lmn_before, :B_lmn_after] .=> unitize(B_unit),
# renamecols=false
# )

function process_events!(events, data, V; kwargs...)
    @chain begin
        process_events!(events, data; kwargs...)
        @rtransform! @astable begin
            V_t = _unitify_V(parent(tview(V, :time)))
            :V = V_t
            # t_us = :t_us_ds[1]
            # t_ds = :t_us_ds[2]
            # :V_us = tview(tview(V, :tstart, t_us), t_us)
            # :V_ds = tview(tview(V, t_ds, :tstop), t_ds)
        end
    end
end


function process_events!(events, data, V, n; kwargs...)
    @chain begin
        process_events!(events, data, V; kwargs...)
        @rtransform! :n = only(tview(n, :time))
    end
end
