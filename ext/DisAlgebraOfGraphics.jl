module DisAlgebraOfGraphics

using AlgebraOfGraphics
using AlgebraOfGraphics.Makie
using AlgebraOfGraphics: density, FigureGrid
import Discontinuity: plot_dist!, plot_dist

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
        p = spec * mapping(m)
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


end