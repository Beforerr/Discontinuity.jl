"""
    HyperbolicTangentFit{T}

A hyperbolic tangent fit.

```math
y = A tanh((x - μ) / σ) + B
```
"""
struct HyperbolicTangentFit{T}
    A::T
    μ::T
    σ::T
    B::T
end

HyperbolicTangentFit(p::AbstractArray{T}) where T = HyperbolicTangentFit{T}(p...)

(fit::HyperbolicTangentFit)(x) = tanh_model(x, fit.A, fit.μ, fit.σ, fit.B)

for func in (:tanh_model, :tanh_model!)
    @eval FitModel(::typeof($func)) = HyperbolicTangentFit
end

function _gradient(fit::HyperbolicTangentFit, x)
    return fit.A / fit.σ * sech((x - fit.μ) / fit.σ)^2
end