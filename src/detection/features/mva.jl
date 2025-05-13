using LsqFit
using LsqFit: rss
using Statistics
using Dates
include("fit.jl")

times(data) = DimensionalData.lookup(dims(data, TimeDim))

function init_p0(data, xdata; σ_min=0)
    x_min, x_max = 0, maximum(xdata)
    A0 = (data[end] - data[1]) / 2
    B0 = (data[end] + data[1]) / 2
    μ0 = xdata[argmin(@. abs(data - B0))]
    # σ0 = (x_max - x_min) / 7
    σ0 = min(x_max - μ0, μ0 - x_min)
    p0 = [A0, μ0, σ0, B0]
    lb = SA_F64[-2abs(A0), x_min, σ_min, B0-abs(A0)] # for `lmfit` keyword argument `lower`, expected AbstractVector{Float64}
    ub = SA_F64[2abs(A0), x_max, Inf, B0+abs(A0)]
    return (p0, lb, ub)
end

"""
    fit_maximum_variance_direction(data, times)

Fit a hyperbolic tangent model to time-series in `data` and extract features
"""
function fit_maximum_variance_direction(data, times; model=tanh_model!, inplace=true)
    # Not enough points
    if length(data) < 4
        return (; t_fit=missing, fit_param=missing, grad=missing, nrmsd=missing, duration=missing)
    end

    dt = Millisecond(1)
    t0 = minimum(times)
    xdata = @. (times - t0) / dt

    p0, lb, ub = init_p0(data, xdata)
    fit = curve_fit(model, xdata, data, p0; lower=lb, upper=ub, inplace)
    p = fit.param

    # Not enough points within 3σ range
    σ_range = (p[2] - 3p[3], p[2] + 3p[3])
    if count(x -> σ_range[1] < x < σ_range[2], xdata) < 3
        return (; t_fit=missing, fit_param=missing, grad=missing, nrmsd=missing, duration=missing)
    end

    μ = p[2]
    f = FitModel(model)(p)
    udt = uconvert(u"s", dt)
    grad = _gradient(f, μ) / udt
    t_fit = t0 + Nanosecond(round(Int, μ * (dt / Nanosecond(1))))
    return (; t_fit, fit_param=p, grad, nrmsd=nrmsd(fit, data), duration=2p[3] * udt)
end

fit_maximum_variance_direction(data) = fit_maximum_variance_direction(parent(data), times(data))

# https://github.com/JuliaNLSolvers/LsqFit.jl/pull/59
# rsquared(data, model, fit) = cor(data, model(t, fit.param))^2


function mva_features(data)
    eigen = mva_eigen(data)
    mva_data = rotate(data, eigen)
    B_l = view(mva_data, :, 1)
    B_n = view(mva_data, :, 3)
    fit = fit_maximum_variance_direction(B_l)
    return (;
        B_lmn_before=parent(mva_data)[1, :],
        B_lmn_after=parent(mva_data)[end, :],
        λ2_over_λ3=eigen.values[2] / eigen.values[3],
        e_max=eigen.vectors[:, 1],
        n_mva=eigen.vectors[:, 3],
        B_n_mva=mean(B_n),
        fit...
    )
end


function mva_transform((V, B), tmin, tmax; tb_min=tmin, tb_max=tmax)
    mva(V(tmin, tmax), B(tb_min, tb_max))
end