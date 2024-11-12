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
        unique!(["t_us", "t_ds"]) # Remove duplicate rows
        select!(Not(cols2remove)) # Remove additional columns
    end
end

"""
    backwards_comp!(df)

Backwards compatibility for old column names
"""
function backwards_comp!(df)
    # renaming
    :"t.d_start" in names(df) && @rename!(df, :t_us = :"t.d_start")
    :"t.d_end" in names(df) && @rename!(df, :t_ds = :"t.d_end")
    :"Vl" in names(df) && @rename!(df, :e_max = :Vl)
    :"V.before" in names(df) && @rename!(df, :V_us = :"V.before")
    :"V.after" in names(df) && @rename!(df, :V_ds = :"V.after")
    :"V.ion.before" in names(df) && @rename!(df, :V_us = :"V.ion.before")
    :"V.ion.after" in names(df) && @rename!(df, :V_ds = :"V.ion.after")
end

"""
    load(path)

Load the data from the given path
"""
load(path) = path |> Arrow.Table |> DataFrame

@kwdef struct DataSet
    name = "events"
    path = missing
    ts = missing
    tau = missing
    detect_func = "detect_variance"
    method = "fit"
end

decode(v::Period) = Dates.format(Time(0) + v, "H:MM:SS")

"""
Return the filename regex of the dataset.
"""
function rfilename(ds::DataSet)
    ts_part = "ts=" * decode(ds.ts)
    tau_part = "tau=" * decode(ds.tau)
    detect_func_part = "detect_func=" * ds.detect_func
    method_part = "method=" * ds.method
    return Regex("^$(ds.name)_tr=.*.*$detect_func_part.*$tau_part.*$ts_part.*$method_part")
end

"""
Return the filename of the dataset using regex.
"""
function filename(ds::DataSet, dir)
    r = rfilename(ds)
    files = [f for f in readdir(dir) if match(r, f) !== nothing]
    if isempty(files)
        error("No file found for $ds")
    elseif length(files) > 1
        error("Multiple files found for $ds: $files")
    end
    return joinpath(dir, files[1])
end


function path(ds::DataSet; dir=".")
    ismissing(ds.path) ? filename(ds, dir) : ds.path
end

load(ds::DataSet; kw...) = (load ∘ path)(ds; kw...)