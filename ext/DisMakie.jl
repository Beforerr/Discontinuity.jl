module DisMakie
using Discontinuity
using Discontinuity: waiting_time
using Beforerr: DistsFit
using Makie
using Distributions
import Discontinuity: plot_fit, plot_fit!, plot_B_mva!
using Dates: AbstractTime, Second
using Discontinuity.SPEDAS: tview, mva, times, tlines!

# log_tickformat = values -> [L"10^{%$(value)}" for value in values]
# axis=(xtickformat = log_tickformat,)

function plot_wt_pdf!(
        fp, τ;
        dists = (Exponential,),
        step = 5,
        xscale = identity,
        yscale = log10,
        add_legend = true,
        legend_title = "Distribution"
    )
    ax = Axis(fp; yscale, xscale, xlabel = "τ (minutes)", ylabel = "p(τ)")
    plot!(ax, DistsFit(τ, dists, step))
    add_legend && axislegend(ax, legend_title; loc = :upperright)
    return fp
end


"""
Plot the waiting time distribution of the data
"""
function plot_wt_pdf(τ::AbstractVector; kwargs...)
    f = Figure()
    plot_wt_pdf!(f[1, 1], τ; kwargs...)
    return f
end

plot_wt_pdf(df; kwargs...) = plot_wt_pdf(waiting_time(df); kwargs...)


function Discontinuity.plot_B_mva!(fp, data, tmin::AbstractTime, tmax::AbstractTime, δt = Second(0); idx = 1:1)
    cs = Makie.wong_colors()

    data_s = tview(data, tmin - δt, tmax + δt)
    B_s = tview(data_s, tmin, tmax)
    B_mva = mva(data_s, B_s)
    time = times(B_mva)
    
    for i in idx
        @views lines!(fp, time, parent(B_mva[:, i]))
    end
    tlines!([tmin, tmax]; linestyle=:dash, color=cs[3])
    
    return fp
end

function Discontinuity.plot_B_mva!(fp, data, event, δt = Second(20))
    return plot_B_mva!(fp, data, event.t_us, event.t_ds, δt)
end

function plot_fit!(fp, data, event, δt)
    plot_B_mva!(fp, data, event, δt)
    time = times(data)
    fit = event.fit
    t_fit = event.t_fit
    t0 = event.t_us
    if !ismissing(fit)
        lines!(fp, time, fit.(time, t0))
        scatter!(fp, t_fit, fit(t_fit, t0); color = :red)
    end
    return fp
end

function plot_fit(data, event, δt = Second(20))
    f = Figure()
    ax = Axis(f[1, 1])
    plot_fit!(ax, data, event, δt)
    return f
end

end
