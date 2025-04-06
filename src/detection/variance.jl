"""
    This module contains functions for detecting discontinuities based on variance analysis.
"""

@kwdef mutable struct VarianceOptions
    std_threshold::Float64 = 2
    fluc_threshold::Float64 = 1
    diff_threshold::Float64 = 0.1
    sparse_threshold::Int = 15
end

norm_std(x; dim=1) = norm(nanstd(x; dim))

"""
    compute_std(data, group_idxs, dim)

Compute standard deviation over a rolling window.
"""
function compute_std(data, group_idxs, ::Val{dim}) where dim
    return tmap(group_idxs) do group_idx
        window_data = selectdim(data, dim, group_idx)
        norm_std(window_data)
    end
end

"""
    compute_combined_std(data, group_idxs, n, ::Val{dim})

Compute combined standard deviation.
"""
function compute_combined_std(data, group_idxs, n, ::Val{dim}) where dim
    len = length(group_idxs)
    return tmap(eachindex(group_idxs)) do i
        if i >= n + 1 && i <= len - n
            group_idx = vcat(group_idxs[i-n], group_idxs[i+n])
            # Slower but memory efficient
            # group_idx = ApplyVector(vcat, group_idxs[i-n], group_idxs[i+n])
            window_data = selectdim(data, dim, group_idx)
            norm_std(window_data)
        else
            missing
        end
    end
end

"""
    compute_index_std!(data::DataFrame, tau; on=:starts)

Compute the standard deviation index based on the given data.

First get the neighbor standard deviations.
"""
function compute_index_std!(data::DataFrame, tau; on=:tstart)
    prev_df = @select(data, :tstart = :tstart .- tau, :std_prev = :std, :len_prev = :len)
    next_df = @select(data, :tstart = :tstart .+ tau, :std_next = :std, :len_next = :len)
    return @chain data begin
        leftjoin!(_, prev_df; on)
        leftjoin!(_, next_df; on)
    end
end

"""Compute the difference index."""
function diff_index(data, ::Val{n}) where n
    T = SVector{n}
    slices = eachrow(data)
    Δ = norm(T(first(slices)) - T(last(slices)))
    m = nanmean(norm.(T.(slices)))
    return Δ / m
end

function diff_index(data, group_idxs, ::Val{dim}) where dim
    n = Val(size(data, 2))
    return tmap(group_idxs) do group_idx
        window_data = selectdim(data, dim, group_idx)
        diff_index(window_data, n)
    end
end

"""
    compute_indices(data, period, n=2; dim=Ti)

Compute all indices based on the given data and tau value.
"""
function compute_indices(data, period, n=2; dim=Ti)
    every = period / n
    times = dims(data, dim).val |> parent
    group_idxs, tstart = groupby_dynamic(times, every, period)
    d = Val(dimnum(data, dim))
    pdata = parent(data)

    len = length.(group_idxs)
    index_diff = diff_index(pdata, group_idxs, d)
    std = compute_std(pdata, group_idxs, d)
    std_combined = compute_combined_std(pdata, group_idxs, n, d)
    df = DataFrame((; tstart, len, std, std_combined, index_diff))

    @chain begin
        compute_index_std!(df, period)
        @rtransform!(
            :time = :tstart + period / 2,
            :tstop = :tstart + period,
            :index_std = :std / max(:std_prev, :std_next),
            :index_fluctuation = :std_combined / (:std_prev + :std_next)
        )
    end
end

"""
    filter_indices!(data::DataFrame;
                   index_std_threshold=2.0,
                   index_fluc_threshold=1.0,
                   index_diff_threshold=0.1,
                   sparse_num=15)

Filter indices to get possible discontinuities.
"""
function filter_indices!(
    data::DataFrame;
    index_std_threshold=2.0,
    index_fluc_threshold=1.0,
    index_diff_threshold=0.1,
    sparse_num=15
)
    @subset! data @byrow begin
        :index_std > index_std_threshold
        isfinite(:index_std)
        :index_fluctuation > index_fluc_threshold
        :index_diff > index_diff_threshold
        :len > sparse_num
        :len_prev > sparse_num
        :len_next > sparse_num
    end
end

"""
    detect_variance(data, tau, n=2; dim=Ti, sparse_num=nothing)

Detect discontinuities based on variance analysis.
"""
function detect_variance(data, tau; n=2, dim=Ti, sparse_num=nothing)
    ts = time_resolution(data; dim)
    sparse_num = isnothing(sparse_num) ? ceil(Int, tau / ts / 3) : sparse_num

    indices = compute_indices(data, tau, n; dim=dim)
    filter_indices!(indices; sparse_num)
end


# https://github.com/brenhinkeller/NaNStatistics.jl/issues/55
# ustrip first provides a performance boost
function compute_indices(data::AbstractArray{Q}, args...; kwargs...) where Q<:Quantity
    return compute_indices(ustrip(data), args...; kwargs...)
end