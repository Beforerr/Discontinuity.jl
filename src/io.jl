using Arrow
using FileIO

DEFAULT_REDUNDANT_COLS = [
    ["B.vec.before.", "B.vec.after."] .* ["l" "m" "n"]
] |> Iterators.flatten

"""
    standardize_df!(df)

Standardize the DataFrame by:
1. Converting Float32 columns to Float64
2. Removing rows with NaN values 
3. Removing duplicate rows based on start/end times
4. Removing redundant columns defined in DEFAULT_REDUNDANT_COLS
"""
function standardize_df!(df)
    cols2remove = names(df) ∩ DEFAULT_REDUNDANT_COLS

    backwards_comp!(df)

    @chain df begin
        transform!(names(df, Float32) .=> ByRow(Float64); renamecols=false) # Convert all columns of Float32 to Float64
        subset!(names(df, Float64) .=> ByRow(isfinite)) # Remove rows with NaN values
        unique!(["t.d_start", "t.d_end"]) # Remove duplicate rows
        select!(Not(cols2remove)) # Remove additional columns
    end
end

"""
    backwards_comp!(df)

Backwards compatibility for old column names
"""
function backwards_comp!(df)
    # renaming
    @rename!(df, :Vl => :e_max, :t_us => "t.d_start", :t_ds => "t.d_end")
end

process!(df::AbstractDataFrame) = begin
    df |>
    dropmissing |>
    keep_good_fit! |>
    standardize_df! |>
    compute_params! |>
    compute_Alfvenicity_params!
end

"""
    load(path)

Load the data from the given path
"""
load(path) = path |> Arrow.Table |> DataFrame

@kwdef struct DataSet
    name = missing
    path = missing
    ts = missing
    tau = missing
    method = "fit"
end

prefix(ds::DataSet) = "updated_events_$(ds.name)_"
suffix(ds::DataSet) = "arrow"

function path(ds::DataSet)
    ismissing(ds.path) ? "$(prefix(ds))$(ds.ts)_$(ds.tau).arrow" : ds.path
end

load(ds::DataSet) = (load ∘ path)(ds)