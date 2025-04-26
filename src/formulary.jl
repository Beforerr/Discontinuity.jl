using Unitful
using Unitful: μ0, Units, mp
using Unitful: BField
using PlasmaFormulary
import PlasmaFormulary: plasma_beta, alfven_velocity, thermal_temperature, NumberDensity

const DEFAULT_B_UNIT = u"nT"
const DEFAULT_L_UNIT = u"km"
const DEFAULT_N_UNIT = u"cm^-3"
const DEFAULT_V_UNIT = u"km/s"
const V_UNIT = u"km/s"
const DEFAULT_T_UNIT = u"eV"

function Alfven_speed(B::BField, n::NumberDensity)
    return B / sqrt(μ0 * n * mp) |> upreferred
end

function gradient_current(dBdt, V)
    return dBdt / (V * μ0) |> upreferred
end

gradient_current(dBdt::Unitful.Frequency, V) =
    gradient_current(dBdt * DEFAULT_B_UNIT, V)

PlasmaFormulary.plasma_beta(T::Real, n::Real, B::Real) = plasma_beta(T * DEFAULT_T_UNIT, n * DEFAULT_N_UNIT, B * DEFAULT_B_UNIT) |> NoUnits
PlasmaFormulary.alfven_velocity(B::Real, n::Real) = alfven_velocity(B * DEFAULT_B_UNIT, n * DEFAULT_N_UNIT) / DEFAULT_V_UNIT |> NoUnits
PlasmaFormulary.thermal_temperature(V::Real, mass=Unitful.mp) = thermal_temperature(V * DEFAULT_V_UNIT, mass)


function Alfven_current(B, n)
    return upreferred(Alfven_speed(B, n) * n * Unitful.q)
end

"""
    safeunitize(x, unit)

Convert the input data `x` to the specified `unit` if `x` does not have the unit.
"""
safeunitize(x, unit) = x
safeunitize(x::Real, unit) = x * unit
safeunitize(x::AbstractArray, unit) = safeunitize.(x, unit)

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

function anisotropy(Bmag::BField, n, T_parp, T_perp)
    (μ0 * n * (T_parp - T_perp) / Bmag^2) |> NoUnits
end

anisotropy(Bmag::BField, n, T3; i=3) = anisotropy(Bmag, n, decompose_T3(T3, i)...)

function decompose_T3(x, i)
    T_parp = x[i]
    T_perp = mean(x[j] for j in eachindex(x) if j != i)
    return T_parp, T_perp
end

anisotropy(B::AbstractVector, args...) = anisotropy(norm(B), args...)

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
