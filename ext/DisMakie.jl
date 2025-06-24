module DisMakie
using Discontinuity: waiting_time
using Beforerr: DistsFit
using Makie
using Distributions
import Discontinuity: plot_fit, plot_fit!
using Dates: Second
using Discontinuity.SPEDAS: tview, mva, times, tlines!

# log_tickformat = values -> [L"10^{%$(value)}" for value in values]
# axis=(xtickformat = log_tickformat,)

function plot_wt_pdf!(fp, τ;
    dists=(Exponential,),
    step=5,
    xscale=identity,
    yscale=log10,
    add_legend=true,
    legend_title="Distribution"
)
    ax = Axis(fp; yscale, xscale, xlabel="τ (minutes)", ylabel="p(τ)")
    plot!(ax, DistsFit(τ, dists, step))
    add_legend && axislegend(ax, legend_title; loc=:upperright)
    return fp
end


"""
Plot the waiting time distribution of the data
"""
function plot_wt_pdf(τ::AbstractVector; kwargs...)
    f = Figure()
    plot_wt_pdf!(f[1, 1], τ; kwargs...)
    f
end

plot_wt_pdf(df; kwargs...) = plot_wt_pdf(waiting_time(df); kwargs...)


function plot_fit!(fp, data, event, δt = Second(20))
    cs = Makie.wong_colors()

    tmin, tmax = event.t_us, event.t_ds
    data_s = tview(data, tmin - δt, tmax + δt)
    B_s = tview(data_s, tmin, tmax)
    B_mva = mva(data_s, B_s)
    time = times(B_mva)
    fit = event.fit
    t_fit = event.t_fit
    t0 = tmin
    B_l = B_mva[:, 1]
    lines!(fp, time, parent(B_l))
    tlines!([tmin, tmax]; linestyle=:dash, color=cs[3])

    if !ismissing(fit)
        lines!(fp, time, fit.(time, t0))
        scatter!(fp, t_fit, fit(t_fit, t0); color=:red)
    end
    fp
end

function plot_fit(data, event, δt = Second(20); )
    f = Figure()
    ax = Axis(f[1, 1])
    plot_fit!(ax, data, event, δt)
    f
end

end