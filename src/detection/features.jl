using SPEDAS: mva, mva_eigen, rotate
using SPEDAS

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
    B_n_cross = sproj(nanmean(data, dims = 1), n_cross)
    return (; n_cross, B_n_cross)
end

function span(x)
    a, b = extrema(x)
    return abs(b - a)
end

function determine_boundary!(df, data; kw...)
    return @rtransform! df @astable begin
        t_us, t_ds = ts_max_distance(tview(data, :tstart, :tstop); kw...)
        :t_us = t_us
        :t_ds = t_ds
    end
end

function process_events!(events, B; kwargs...)
    return @chain events begin
        remove_duplicates(_)
        @rtransform! @astable begin
            :B_us = SV3(tselect(B, :t_us))
            :B_ds = SV3(tselect(B, :t_ds))
            :ω = angle_between(:B_us, :B_ds)
        end
        @rtransform! $AsTable = cross_features(parent(tview(B, :t_us, :t_ds)))
        @rtransform! $AsTable = stat_features(tview(B, :t_us, :t_ds))
        @rtransform! $AsTable = mva_features(tview(B, :t_us, :t_ds))
    end
end

function process_events_V!(events, V, δt; verbose = false)
    return @rtransform! events @astable begin
        t_us = :t_us
        t_ds = :t_ds
        :V = SV3(tselect(V, :time))

        Vs_us = tview(V, t_us - δt, t_us)
        Vs_ds = tview(V, t_ds, t_ds + δt)
        flag = length(Vs_us) == 0 || length(Vs_ds) == 0

        # calculate the largest change in the l direction
        :ΔV_l_max = flag ? missing : begin
                t_V_us = tselect(times(Vs_us), t_us)
                t_V_ds = tselect(times(Vs_ds), t_ds)
                V_subset = tview(V, t_V_us, t_V_ds)
                V_subset_l = sproj.(eachslice(V_subset; dims = TimeDim), Ref(:e_max))
                span(V_subset_l)
            end

        verbose && length(Vs_us) == 0 && @info "No upstream data between $(t_us - δt) and $t_us"
        verbose && length(Vs_ds) == 0 && @info "No downstream data between $t_ds and $(t_ds + δt)"
        :V_us = length(Vs_us) == 0 ? missing : SV3(tselect(Vs_us, t_us))
        :V_ds = length(Vs_ds) == 0 ? missing : SV3(tselect(Vs_ds, t_ds))
    end
end

function process_events!(events, B, V; δt = Minute(1), kwargs...)
    return @chain events begin
        process_events!(_, B; kwargs...)
        process_events_V!(_, V, δt)
    end
end


function process_events!(events, B, V, n; kwargs...)
    return @chain begin
        process_events!(events, B, V; kwargs...)
        @rtransform! :n = only(tselect(n, :time))
    end
end
