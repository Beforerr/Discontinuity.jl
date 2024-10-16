using Unitful
using Unitful: μ0, Units
using PlasmaFormulary

DEFAULT_B_UNIT = u"nT"
DEFAULT_L_UNIT = u"km"
DEFAULT_N_UNIT = u"cm^-3"
DEFAULT_V_UNIT = u"km/s"
DEFAULT_T_UNIT = u"eV"

"""
    unitize([f], unit::Units)

Convert the input data `f` to the specified `unit` if the data type is `Float64`.
"""
function unitize(f, unit::Units)
    return eltype(f) <: Float64 ? f .* unit : f
end

unitize(unit::Units) = f -> unitize(f, unit)

function unitize!(df, col, unit)
    # check if the type of the column data is float
    df[!, col] = unitize(df[:, col], unit)
end

function unitize!(
    df;
    Bcols=[DEFAULT_B_COL],
    ncols=[DEFAULT_N_COL],
    Tcols=[DEFAULT_T_COL],
    Lcols=DEFAULT_L_COLS,
    B_unit=DEFAULT_B_UNIT,
    n_unit=DEFAULT_N_UNIT,
    T_unit=DEFAULT_T_UNIT,
    L_unit=DEFAULT_L_UNIT,
)
    transform!(
        df,
        Bcols .=> unitize(B_unit),
        ncols .=> unitize(n_unit),
        Tcols .=> unitize(T_unit),
        Lcols .=> unitize(L_unit);
        renamecols=false
    )
end


function anisotropy(B, density, para_temp, perp_temp)
    Λ = @. (μ0 * density * (para_temp - perp_temp) / B^2) |> NoUnits
    return Λ
end

function calc_beta!(
    df;
    B=DEFAULT_B_COL,
    n=DEFAULT_N_COL,
    T=DEFAULT_T_COL,
)
    T = df[:, T]
    n = df[:, n]
    B = df[:, B]
    df[!, :β] = PlasmaFormulary.beta.(T, n, B)
    return df
end


function calc_T!(df, V; kwargs...)
    V = df[:, V]
    mass = Unitful.mp
    df[!, :T] = PlasmaFormulary.thermal_temperature.(V, mass; kwargs...)
    return df
end


"""
Cos rotation angle between two vectors.
"""
function cos_rotation_angle(v1, v2)
    return v1 ⋅ v2 / (norm(v1) * norm(v2))
end

"""Rotation angle between two vectors in radians."""
rotation_angle(v1, v2) = acos(cos_rotation_angle(v1, v2))
"""Rotation angle between two vectors in degrees."""
rotation_angle_d(v1, v2) = rad2deg(rotation_angle(v1, v2))