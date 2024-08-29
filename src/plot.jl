log_tickformat = values -> [L"10^{%$(value)}" for value in values]
# axis=(xtickformat = log_tickformat,)
using AlgebraOfGraphics
using AlgebraOfGraphics: density

"""
Plot the density distribution of the thickness and current density (default) of the data layer.
"""
function plot_dist(
    data_layer;
    maps=[l_log_map, l_norm_log_map, j_log_map, j_norm_log_map],
    axis=(yscale=log10,),
    datalimits=extrema,
    figure_kwargs=(size=(1200, 300),),
    visual=visual(Lines)
)

    plt = data_layer * density(datalimits=datalimits)
    plt = plt * visual
    plts = [plt * mapping(m) for m in maps]

    fig = Figure(; figure_kwargs...)
    axs = [fig[i, j] for j in 1:size(maps, 2) for i in 1:size(maps, 1)]
    grids = map(axs, plts) do ax, p
        draw!(ax, p, axis=axis)
    end
    add_labels!(axs)
    pretty_legend!(fig, grids[1])

    fig
end