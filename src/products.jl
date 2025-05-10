using SPEDAS.SpaceDataModel: AbstractProduct
using SPEDAS: rotate

struct FitProduct{F,P} <: AbstractProduct
    f::F
    param::P
    t0::DateTime
end

function (p::FitProduct)(t0, t1)
    tspan = DateTime(t0):Millisecond(10):DateTime(t1)
    xspan = @. (tspan - p.t0) / Millisecond(1)
    data = p.f.(xspan, Ref(p.param))
    return DimVector(data, Ti(tspan))
end

struct MVAProduct{T,R} <: AbstractProduct
    V::T
    rotation::R
end

MVAProduct(V, B, tspan) = MVAProduct(V, mva_eigen(B(tspan...)))
MVAProduct(V) = MVAProduct(V, V)

function (p::MVAProduct)(t0, t1)
    data = p.V(t0, t1)
    p.rotation isa Eigen ? rotate(data, p.rotation) : mva(data, p.rotation(t0, t1))
end