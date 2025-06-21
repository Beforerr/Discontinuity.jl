using Unitful
using Unitful: μ0, Units, mp, q
using Unitful: BField, Velocity, Temperature, Frequency, Energy
import PlasmaFormulary
import PlasmaFormulary: plasma_beta, NumberDensity, CurrentDensity

const qe = Unitful.q
const B_UNIT = u"nT"
const L_UNIT = u"km"
const N_UNIT = u"cm^-3"
const J_UNIT = u"nA/m^2"
const V_UNIT = u"km/s"
const T_UNIT = u"eV"
const TE_UNIT = u"eV"
const QuantityLikeType = Union{Quantity,AbstractArray{<:Quantity}}

_unitify_L(x) = safeunitize(x, L_UNIT)
_unitify_B(x) = safeunitize(x, B_UNIT)
_unitify_n(x) = safeunitize(x, N_UNIT)
_unitify_V(x) = safeunitize(x, V_UNIT)
_unitify_T(x) = safeunitize(x, T_UNIT)

uless(x) = NoUnits(x)
uless(x::NumberDensity) = NoUnits(x / N_UNIT)
uless(x::CurrentDensity) = NoUnits(x / J_UNIT)
uless(x::Unitful.Length) = NoUnits(x / L_UNIT)
uless(x::Velocity) = NoUnits(x / V_UNIT)
uless(x::Temperature) = NoUnits(x * Unitful.k / TE_UNIT)
uless(x::Energy) = NoUnits(x / TE_UNIT)

function gradient_current(dBdt, V)
    return dBdt / (V * μ0) |> upreferred
end

gradient_current(dBdt::Frequency, V) =
    gradient_current(dBdt * B_UNIT, V)

for f in (:Alfven_velocity, )
    @eval $f(args::Vararg{QuantityLikeType}) = PlasmaFormulary.$f(args...)
end

Alfven_velocity(B, n) = Alfven_velocity(_unitify_B(B), _unitify_n(n))
plasma_beta(T, n, B) = plasma_beta(_unitify_T(T), _unitify_n(n), _unitify_B(B))

inertial_length(n, q=qe, m=mp) = PlasmaFormulary.inertial_length(_unitify_n(n), q, m)
thermal_temperature(V, mass=Unitful.mp) = PlasmaFormulary.thermal_temperature(_unitify_V(V), mass)
Alfven_speed(B, n) = PlasmaFormulary.Alfven_speed(_unitify_B(B), _unitify_n(n))
Alfven_current(B, n) = upreferred(Alfven_speed(B, n) * _unitify_n(n) * q)

"""
    safeunitize(x, unit)

Convert the input data `x` to the specified `unit` if `x` does not have the unit.
"""
safeunitize(x, unit) = eltype(x) <: Real ? x * unit : x

unitize(unit::Units) = f -> safeunitize(f, unit)

unitize!(df, col, unit) = (df[!, col] = safeunitize(df[!, col], unit))

function unitize!(
    df;
    Bcols=[DEFAULT_B_COL],
    ncols=[DEFAULT_N_COL],
    Vcols=[],
    Tcols=[DEFAULT_T_COL],
    Lcols=DEFAULT_L_COLS ∩ names(df),
)
    transform!(
        df,
        Bcols .=> _unitify_B,
        ncols .=> _unitify_n,
        Vcols .=> _unitify_V,
        Tcols .=> _unitify_T,
        Lcols .=> _unitify_L;
        renamecols=false
    )
end

function anisotropy(Bmag::BField, n, T_parp, T_perp)
    (μ0 * n * (T_parp - T_perp) / Bmag^2) |> NoUnits
end

anisotropy(Bmag, n, T3; i=3) = anisotropy(Bmag, n, decompose_T3(T3, i)...)

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

# Unitize the input arguments

anisotropy(Bmag, n , T_parp, T_perp) = anisotropy(
    _unitify_B(Bmag),
    _unitify_n(n),
    _unitify_T(T_parp),
    _unitify_T(T_perp)
)
