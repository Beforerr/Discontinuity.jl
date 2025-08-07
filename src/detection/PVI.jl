function PVI_deriv(B; dim = nothing)
    dim = @something dim dimnum(B, TimeDim)
    # Compute increments
    odim = other_dims(B, dim)
    out = sum(abs2 ∘ ustrip, tderiv(TimeseriesUtilities.DiffQ, B; dim); dims = odim)
    rms = sqrt(mean(out))
    return out .= sqrt.(out) ./ rms
end

function _PVI(B, n, ::Val{dim}) where {dim}
    # Compute increments
    N = size(B, dim)
    ΔB = selectdim(B, dim, (n + 1):N) - selectdim(B, dim, 1:(N - n))
    # Magnitudes of vector increments
    odim = other_dims(B, dim)
    out = sum(abs2, ΔB; dims = odim)
    rms = sqrt(mean(out))
    out .= sqrt.(out) ./ rms
    return vec(out)
end


# TODO: add time interpolation

"""
    PVI(B, n, ::Val{dim})

Compute the Partial Variance of Increments (PVI) for a time series `B` along dimension `dim` with lag `n`.
Standard deviations σi of magnetic field increments are computed `every` time units.

```math
\\mathrm{PVI}_{t, τ} = \\frac{|Δ \\mathbf{B}(t, τ)|}{\\sqrt{\\langle | Δ \\mathbf{B}(t, τ) |^2\\rangle}}
```

# Returns a Vector of length N - τ containing PVI values.
"""
function PVI(B::AbstractDimArray, τ = time_resolution(B); every = nothing, dim = nothing)
    dim = @something dim dimnum(B, TimeDim)
    metadata = (; label = "PVI", unit = "")
    n = ceil(Int, τ / time_resolution(B))
    @info "PVI lag: $n"
    return if isnothing(every)
        data = _PVI(parent(B), n, Val(dim))
        dims = (DimensionalData.dims(B, dim)[1:(end - n)],)
        rebuild(B; data, dims, metadata)
    else
        group_idxs, idxs = groupby_dynamic(times(B), every)
        mapreduce(vcat, group_idxs) do group_idx
            B_s = selectdim(B, dim, group_idx)
            data = _PVI(parent(B_s), n, Val(dim))
            dims = (DimensionalData.dims(B_s, dim)[1:(end - n)],)
            rebuild(B_s; data, dims, metadata)
        end
    end
end


"""
    group_seeds_into_intervals(seed_indices, gap)

Group an ordered vector of seed indices into contiguous intervals.
Seeds separated by no more than `gap` belong to the same interval.
Returns a vector of (start_idx, end_idx) tuples.
"""
function group_seeds_into_intervals(seed_indices, gap)
    intervals = Tuple{Int, Int}[]
    if isempty(seed_indices)
        return intervals
    end
    start_idx = seed_indices[1]
    prev_idx = start_idx
    for idx in seed_indices[2:end]
        if idx - prev_idx <= gap
            prev_idx = idx
        else
            push!(intervals, (start_idx, prev_idx))
            start_idx = idx
            prev_idx = idx
        end
    end
    push!(intervals, (start_idx, prev_idx))
    return intervals
end

function coherent_intervals(pvi, times; seed_threshold = 4.0, gap = 3)
    # Step 1: Identify seeds
    seed_indices = findall(pvi .>= seed_threshold)

    # Step 2: Group seeds into intervals
    intervals = group_seeds_into_intervals(seed_indices, gap)

    # Step 3: Expand intervals to local minima below boundary_threshold
    coherent_intervals = Tuple{Int, Int}[]
    for (start_idx, end_idx) in intervals
        # Find local minima on both sides using adaptive threshold
        start_threshold = pvi[start_idx] / 8
        end_threshold = pvi[end_idx] / 8
        interval_start = find_local_minimum(pvi, start_idx, -1, start_threshold)
        interval_end = find_local_minimum(pvi, end_idx, 1, end_threshold)
        push!(coherent_intervals, (interval_start, interval_end))
    end
    return map(unique(coherent_intervals)) do (i0, i1)
        # Find time with highest PVI index within the interval
        max_pvi_idx = @views argmax(pvi[i0:i1]) + i0 - 1
        (tstart = times[i0], time = times[max_pvi_idx], tstop = times[i1])
    end
end

coherent_intervals(pvi::AbstractDimArray; kw...) = coherent_intervals(parent(pvi), times(pvi); kw...)

mul(Δt, ratio) = Millisecond(round(Int, Dates.value(Δt) * ratio))

function compute_index_std(x, tmin, tmax; Δt_ratio = 1)
    Δt = tmax - tmin
    std = norm_std(parent(tview(x, tmin, tmax)))
    left_std = norm_std(tview(x, tmin - mul(Δt, Δt_ratio), tmin))
    right_std = norm_std(tview(x, tmax, tmax + mul(Δt, Δt_ratio)))
    return std / max(left_std, right_std)
end

@inline relative_variation(x; dim = 1) =
    norm_std(x; dim) / mean(norm.(eachslice(x; dims = dim)))

"""
Relative variation σ/⟨B⟩ at left and right boundaries
"""
function compute_index_stability(x, tmin, tmax; Δt_ratio = 1)
    Δt = tmax - tmin
    x_up = tview(x, tmin - mul(Δt, Δt_ratio), tmin)
    x_down = tview(x, tmax, tmax + mul(Δt, Δt_ratio))
    return max(relative_variation.((x_down, x_up))...)
end

function compute_index_stability!(df, x; stability_threshold = 0.3, kw...)
    initial_count = size(df, 1)
    @chain df begin
        @rtransform! :index_stability = compute_index_stability(x, :tstart, :tstop; kw...)
        @subset! :index_stability .<= stability_threshold
    end
    stability_count = size(df, 1)
    @info "After stability filter: $(stability_count) intervals (filtered out: $(initial_count - stability_count))"
    return df
end

# filter: remove wave-like structures
function compute_index_std!(df, x; tmin = :t_us, tmax = :t_ds, std_threshold = 2, kw...)
    initial_count = size(df, 1)
    @chain df begin
        @rtransform! :index_std = compute_index_std(x, $tmin, $tmax; kw...)
        @subset! :index_std .> std_threshold
    end
    std_count = size(df, 1)
    @info "After std filter: $(std_count) intervals (filtered out: $(initial_count - std_count))"
    return df
end

function compute_index_consistency(x, tmin, tmax; N = 2, kw...)
    duration = tmax - tmin
    x_left = tselect(x, tmin)
    x_right = tselect(x, tmax)
    
    # Generate N equally spaced points between tmin and tmax
    interior_angles = map(1:N) do i
        t_i = tmin + mul(duration, i / (N + 1))
        x_i = tselect(x, t_i)
        max(angle_between(x_i, x_left), angle_between(x_i, x_right))
    end
    
    return angle_between(x_left, x_right) / maximum(interior_angles)
end

# filter: remove hole-like and wave-like structures
function compute_index_consistency!(df, x; tmin = :t_us, tmax = :t_ds, consistency_threshold = 1, kw...)
    initial_count = size(df, 1)
    @chain df begin
        @rtransform! :index_consistency = compute_index_consistency(x, $tmin, $tmax; kw...)
        @subset! :index_consistency .>= consistency_threshold
    end
    consistency_count = size(df, 1)
    @info "After consistency filter: $(consistency_count) intervals (filtered out: $(initial_count - consistency_count))"
    return df
end

"""
    select_current_sheets(intervals, data; min_width_factor=0.5, stability_threshold=0.3)

Select intervals that qualify as current sheets (CS) based on:
- Relatively stable magnetic field in the boundaries (not shorter than the temporal half-thickness)
- Reverse sign: angle between two boundaries is larger than angle between boundaries and the center of the interval

# Arguments
- `intervals`: Vector of (tmin, tmax) time intervals
- `std_threshold`: Maximum relative variation allowed in boundary regions
"""
function select_current_sheets!(df, x; std_threshold = 2, stability_threshold = 0.3, consistency_threshold = 1)
    initial_count = size(df, 1)
    compute_index_std!(df, x; std_threshold)
    compute_index_stability!(df, x; stability_threshold)
    compute_index_consistency!(df, x; consistency_threshold)
    final_count = size(df, 1)
    @info "Total CS candidates: $(final_count) out of $(initial_count) initial intervals"
    return df
end
