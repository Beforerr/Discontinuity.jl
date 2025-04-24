using LsqFit
using Statistics
using Dates
using Enzyme
using Enzyme: Reverse, gradient
include("fit.jl")

times(data) = DimensionalData.lookup(dims(data, TimeDim))


"""
    fit_maximum_variance_direction(data, times)

Fit a hyperbolic tangent model to time-series in `data` and extract features
"""
function fit_maximum_variance_direction(data, times; model=tanh_model!, inplace=true)
    # Not enough points
    if length(data) < 4
        return (t_fit=missing, fit_param=missing)
    end

    t0 = minimum(times)
    tspan = maximum(times) - t0
    xdata = @. (times - t0) / tspan
    x_min, x_max = 0, 1

    A0 = (data[end] - data[1]) / 2
    B0 = (data[end] + data[1]) / 2
    μ0 = x_min + (x_max - x_min) / 2
    σ0 = (x_max - x_min) / 7

    p0 = [A0, μ0, σ0, B0]
    lb = SA[-2abs(A0), x_min, 0.0, B0-abs(A0)]
    ub = SA[2abs(A0), x_max, Inf, B0+abs(A0)]
    fit = curve_fit(model, xdata, data, p0; lower=lb, upper=ub, inplace)
    p = fit.param

    μ = p[2]
    f = FitModel(model)(p)
    grad = gradient(Reverse, f, μ)[1] / uconvert(u"s", tspan)
    t_fit = t0 + Nanosecond(round(Int, μ * (tspan / Nanosecond(1))))
    return (; t_fit, fit_param=p, grad)
end

fit_maximum_variance_direction(data) = fit_maximum_variance_direction(parent(data), times(data))

# https://github.com/JuliaNLSolvers/LsqFit.jl/pull/59
# rsquared(data, model, fit) = cor(data, model(t, fit.param))^2


function mva_features(data)
    eigen = mva_eigen(data)
    mva_data = rotate(data, eigen)
    B_l = view(mva_data, :, 1)
    fit = fit_maximum_variance_direction(B_l)
    return (;
        B_lmn_before=parent(mva_data)[1, :],
        B_lmn_after=parent(mva_data)[end, :],
        λ2_over_λ3=eigen.values[2] / eigen.values[3],
        e_max=eigen.vectors[:, 1],
        n_mva=eigen.vectors[:, 3],
        B_n=@views mean(mva_data[:, 3]),
        fit...
    )
end