log_tickformat = values -> [L"10^{%$(value)}" for value in values]
# axis=(xtickformat = log_tickformat,)
using AlgebraOfGraphics
using AlgebraOfGraphics: density, FigureGrid
using Beforerr: DistsFit

FIGURE_KWARGS = (size=(1200, 300),)

"""
Plot the density distribution of the layer.
"""
function plot_dist!(
    axs,
    l::Layer,
    maps;
    axis=(yscale=log10,),
    datalimits=extrema,
    visual=visual(Lines),
    draw_kwargs=Dict()
)
    spec = l * density(datalimits=datalimits) * visual
    return map(axs, maps) do ax, m
        p = spec * mapping(Pair(m...))
        draw!(ax, p; axis=axis, draw_kwargs...)
    end
end


function plot_dist!(
    f::Figure,
    l::Layer,
    maps;
    plotopts=PlotOpts(),
    kwargs...
)
    axs = [f[i, j] for i in 1:size(maps, 1) for j in 1:size(maps, 2)]
    grids = plot_dist!(axs, l, maps; kwargs...)
    fg = FigureGrid(f, grids[1])
    process_opts!(fg, axs, plotopts)
    return f
end

plot_dist!(l::Layer, maps; kwargs...) = plot_dist!(current_figure(), l, maps; kwargs...)

"""Plot the density distribution of the layer."""
plot_dist(l::Layer, maps; figure=(;), kwargs...) = plot_dist!(Figure(; figure...), l, maps; kwargs...)

plot_dist(l::Layer; maps=[l_log_map l_norm_log_map j_log_map j_norm_log_map], kwargs...) = plot_dist(l, maps; figure=FIGURE_KWARGS, kwargs...)
@deprecate plot_dist(l::Layer; maps, kw...) plot_dist(l::Layer, maps; kw...)

function waiting_time(time; δt=Dates.Minute(1))
    # unique and order the time
    τ = diff(time |> unique |> sort)
    return τ ./ δt
end

waiting_time(df::AbstractDataFrame, col=:time; kwargs...) = waiting_time(df[:, col])

"""
Plot the waiting time distribution of the data
"""
function plot_wt_pdf(
    τ;
    dists=(Exponential,),
    step=5,
    xscale=identity
)
    f = Figure()
    ax = Axis(f[1, 1], yscale=log10, xscale=xscale, xlabel="τ (minutes)", ylabel="p(τ)")
    plot!(ax, DistsFit(τ, dists, step))
    # add the legend (fit parameters)
    axislegend(ax, loc=:upperright)
    return f
end

plot_wt_pdf(df::AbstractDataFrame; kwargs...) = plot_wt_pdf(waiting_time(df); kwargs...)