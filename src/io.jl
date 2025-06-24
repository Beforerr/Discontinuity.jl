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
        transform!(names(df, Union{Float32,Missing}) .=> ByRow(Float64); renamecols=false) # Convert all columns of Float32 to Float64
        subset!(names(df, Union{Float64,Missing}) .=> ByRow(isfinite); skipmissing=true) # Remove rows with NaN values
        remove_duplicates()
        select!(Not(cols2remove)) # Remove additional columns
    end
end

"""
    remove_duplicates(df)

Remove duplicates in the DataFrame.

Mark for deletion if either:
1. Current row has same start/end time as previous row
2. Current event starts before previous event ends (with this, checking 1 is redundant) 
"""
@views function remove_duplicates(df)
    sort!(df, :t_us) # Ensure proper time ordering
    mask = trues(nrow(df)) # Preallocate
    t_us = df.t_us
    t_ds = df.t_ds
    for i in 2:nrow(df) # Start from second row, compare with previous row
        mask[i] = t_us[i] > t_ds[i-1]
    end
    return df[mask, :]
end


"""
    backwards_comp!(df)

Backwards compatibility for old column names
"""
function backwards_comp!(df)
    cols = names(df)
    # renaming
    :"t.d_start" in cols && @rename!(df, :t_us = :"t.d_start")
    :"t.d_end" in cols && @rename!(df, :t_ds = :"t.d_end")
    :"Vl" in cols && @rename!(df, :e_max = :Vl)
    :"Vn" in cols && @rename!(df, :e_min = :Vn)
    "e_min" in cols && @rename!(df, :e_min = :n_mva)
    :"V.before" in cols && @rename!(df, :V_us = :"V.before")
    :"V.after" in cols && @rename!(df, :V_ds = :"V.after")
    :"V.ion.before" in cols && @rename!(df, :V_us = :"V.ion.before")
    :"V.ion.after" in cols && @rename!(df, :V_ds = :"V.ion.after")
end

"""
    load(path)

Load the data from the given path
"""
load(path) = path |> Arrow.Table |> DataFrame

@kwdef struct DataSet
    name = "events"
    path = nothing
    ts = nothing
    tau = nothing
    detect_func = "detect_variance"
    method = "fit"
end

decode(v::Period) = Dates.format(Time(0) + v, "H:MM:SS")

"""
Return the filename regex of the dataset.
"""
function rfilename(ds::DataSet)
    ts_part = isnothing(ds.ts) ? "" : "ts=" * decode(ds.ts)
    tau_part = isnothing(ds.tau) ? "" : "tau=" * decode(ds.tau)
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
        error("No file found for $ds with regex pattern: $r, in $dir")
    elseif length(files) > 1
        error("Multiple files found for $ds: $files")
    end
    return joinpath(dir, files[1])
end


path(ds::DataSet; dir=".") = something(ds.path, filename(ds, dir))

load(ds::DataSet; kw...) = (load ∘ path)(ds; kw...)