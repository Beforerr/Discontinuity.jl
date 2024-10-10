using Arrow
using FileIO

function standardize_df!(df)
    @chain df begin
        transform!(names(df, Float32) .=> ByRow(Float64); renamecols=false) # Convert all columns of Float32 to Float64
        subset!(names(df, Float64) .=> ByRow(isfinite)) # Remove rows with NaN values
        unique!(["t.d_start", "t.d_end"]) # Remove duplicate rows
    end
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