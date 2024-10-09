
using FileIO

function standardize_df!(df)
    @chain df begin
        transform!(names(df, Float32) .=> ByRow(Float64); renamecols=false) # Convert all columns of Float32 to Float64
        subset!(names(df, Float64) .=> ByRow(isfinite)) # Remove rows with NaN values
        unique!(["t.d_start", "t.d_end"]) # Remove duplicate rows
    end
end

process!(df::AbstractDataFrame) = df |> keep_good_fit! |> standardize_df! |> compute_params! |> compute_Alfvenicity_params!

"""
    load(path)

Load the data from the given path and process it.
"""
function load(path)
    df = path |> Arrow.Table |> DataFrame |> dropmissing
    df |> process!
end

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

function load(ds::DataSet)
    df = (load ∘ path)(ds)
end