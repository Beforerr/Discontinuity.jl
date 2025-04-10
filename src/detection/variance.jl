"""
    This module contains functions for detecting discontinuities based on variance analysis.
"""

@kwdef mutable struct VarianceOptions
    std_threshold::Float64 = 2
    fluc_threshold::Float64 = 1
    diff_threshold::Float64 = 0.1
    sparse_threshold::Int = 15
end

norm_std(x; dims=1) = norm(nanstd(x; dims))

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
function compute_combined_std(data, idx1::UnitRange, idx2::UnitRange, ::Val{dim}) where dim
    # Slower but memory efficient
    # group_idx = ApplyVector(vcat, idx1, idx2)
    group_idx = vcat(idx1, idx2)
    window_data = selectdim(data, dim, group_idx)
    norm_std(window_data)
end

compute_combined_std(data, idx1s, idx2s, d) =
    tmap(eachindex(idx1s)) do i
        compute_combined_std(data, idx1s[i], idx2s[i], d)
    end

function compute_index_fluctuation!(df, data, d; fluc_threshold=1)
    @chain df begin
        @transform! :std_combined = compute_combined_std(data, :group_idx_prev, :group_idx_next, d)
        @transform! :index_fluctuation = @. :std_combined / (:std_prev + :std_next)
        @subset! :index_fluctuation .> fluc_threshold
    end
end

"""
    compute_index_std!(df; std_threshold=2)

Compute the standard deviation index based on the given data.

First get the neighbor standard deviations.
"""
function compute_index_std!(df; std_threshold=2)
    return @chain df begin
        @transform! :index_std = @. :std / max(:std_prev, :std_next)
        @subset! :index_std .> std_threshold
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

function compute_index_diff!(df, data, d; diff_threshold=0.1)
    @chain df begin
        @transform!(:index_diff = diff_index(data, :group_idx, d))
        @subset!(:index_diff .> diff_threshold)
    end
end

"""
    compute_indices(data, period, n=2; dim=Ti)

Compute all indices based on the given data and tau value.
"""
function compute_indices(data, period, sparse_num; n=2, dim=Ti, std_threshold=2, fluc_threshold=1, diff_threshold=0.1)
    every = period / n
    times = dims(data, dim).val |> parent
    d = Val(dimnum(data, dim))
    pdata = parent(data)

    group_idx, tstart = groupby_dynamic(times, every, period)
    len = length.(group_idx)
    std = compute_std(pdata, group_idx, d)

    df = DataFrame((; tstart, len, std, group_idx))
    on = :tstart
    prev_df = @select(df, $on = $on .- period, :std_prev = :std, :len_prev = :len, :group_idx_prev = :group_idx)
    next_df = @select(df, $on = $on .+ period, :std_next = :std, :len_next = :len, :group_idx_next = :group_idx)

    @chain df begin
        leftjoin!(_, prev_df; on)
        leftjoin!(_, next_df; on)
        dropmissing!
        @subset! begin
            :len .> sparse_num
            :len_prev .> sparse_num
            :len_next .> sparse_num
        end
        compute_index_std!(_; std_threshold)
        compute_index_diff!(_, pdata, d; diff_threshold)
        compute_index_fluctuation!(_, pdata, d; fluc_threshold)
        @transform!(
            :time = :tstart .+ period / 2,
            :tstop = :tstart .+ period,
        )
    end
end

"""
    detect_variance(data, tau, n=2; dim=Ti, sparse_num=nothing)

Detect discontinuities based on variance analysis.
"""
function detect_variance(data, tau; n=2, dim=Ti, sparse_num=nothing, kwargs...)
    ts = time_resolution(data)
    sparse_num = isnothing(sparse_num) ? ceil(Int, tau / ts / 3) : sparse_num
    compute_indices(data, tau, sparse_num; n, dim, kwargs...)
end


# https://github.com/brenhinkeller/NaNStatistics.jl/issues/55
# ustrip first provides a performance boost
function compute_indices(data::AbstractArray{Q}, args...; kwargs...) where Q<:Quantity
    return compute_indices(ustrip(data), args...; kwargs...)
end