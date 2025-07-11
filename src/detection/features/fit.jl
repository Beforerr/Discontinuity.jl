"""Auto differentiation is used when `Enzyme` is available."""
function _gradient end

@. tanh_model(x, A, μ, σ, B) = A * tanh((x - μ) / σ) + B
tanh_model!(y, x, A, μ, σ, B) = y .= tanh_model.(x, A, μ, σ, B)

function jacobian_tanh!(J, x, A, μ, σ, B)
    @. J[:, 1] = tanh((x - μ) / σ)           # ∂f/∂A
    @views @. J[:, 2] = -A / σ * (1 - J[:, 1]^2)  # ∂f/∂μ
    @views @. J[:, 3] = -A * (x - μ) / σ^2 * (1 - J[:, 1]^2) # ∂f/∂σ
    @. J[:, 4] = 1.0                               # ∂f/∂B
end
jacobian_tanh!(J, x, p) = jacobian_tanh!(J, x, p...)

tanh_model!(y, x, p) = tanh_model!(y, x, p...)
tanh_model(x, p) = tanh_model(x, p...)

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

HyperbolicTangentFit(p) = HyperbolicTangentFit(p...)

(fit::HyperbolicTangentFit)(x) = tanh_model(x, fit.A, fit.μ, fit.σ, fit.B)
(fit::HyperbolicTangentFit)(t, t0) = fit((t - t0) / Millisecond(1))

for func in (:tanh_model, :tanh_model!)
    @eval FitModel(::typeof($func)) = HyperbolicTangentFit
end

function _gradient(fit::HyperbolicTangentFit, x)
    return fit.A / fit.σ * sech((x - fit.μ) / fit.σ)^2
end


# normalized root mean squared deviation (NRMSD)
function nrmsd(fit, data)
    rmse = sqrt(rss(fit) / length(data))
    amin, amax = extrema(data)
    return rmse / (amax - amin)
end