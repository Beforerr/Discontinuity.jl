using Unitful
using Unitful: μ0, Units
using PlasmaFormulary
import PlasmaFormulary: plasma_beta, alfven_velocity, thermal_temperature

const DEFAULT_B_UNIT = u"nT"
const DEFAULT_L_UNIT = u"km"
const DEFAULT_N_UNIT = u"cm^-3"
const DEFAULT_V_UNIT = u"km/s"
const DEFAULT_T_UNIT = u"eV"

PlasmaFormulary.plasma_beta(T::Real, n::Real, B::Real) = plasma_beta(T * DEFAULT_T_UNIT, n * DEFAULT_N_UNIT, B * DEFAULT_B_UNIT) |> NoUnits
PlasmaFormulary.alfven_velocity(B::Real, n::Real) = alfven_velocity(B * DEFAULT_B_UNIT, n * DEFAULT_N_UNIT) / DEFAULT_V_UNIT |> NoUnits
PlasmaFormulary.thermal_temperature(V::Real, mass=Unitful.mp) = thermal_temperature(V * DEFAULT_V_UNIT, mass)

"""
    unitize([f], unit::Units)

Convert the input data `f` to the specified `unit` if the data type is `Float64`.
"""
function safeunitize(f, unit::Units)
    return eltype(f) <: Union{Missing,Float64} ? f .* unit : f
end

unitize(unit::Units) = f -> safeunitize(f, unit)

function unitize!(df, col, unit)
    # check if the type of the column data is float
    df[!, col] = safeunitize(df[:, col], unit)
end

function unitize!(
    df;
    Bcols=[DEFAULT_B_COL],
    ncols=[DEFAULT_N_COL],
    Vcols=[],
    Tcols=[DEFAULT_T_COL],
    Lcols=DEFAULT_L_COLS ∩ names(df),
    B_unit=DEFAULT_B_UNIT,
    n_unit=DEFAULT_N_UNIT,
    V_unit=DEFAULT_V_UNIT,
    T_unit=DEFAULT_T_UNIT,
    L_unit=DEFAULT_L_UNIT,
)
    transform!(
        df,
        Bcols .=> unitize(B_unit),
        ncols .=> unitize(n_unit),
        Vcols .=> unitize(V_unit),
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
    @transform!(df, :β = passmissing(PlasmaFormulary.plasma_beta).($T, $n, $B))
end


function calc_T!(df, V; kwargs...)
    return @transform!(df, :T = passmissing(thermal_temperature).($V; kwargs...))
end
