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
    figure_kwargs=(size=(1200, 300),)
)

    plt = data_layer * density(datalimits=datalimits) * visual(Lines)
    plts = [plt * mapping(m) for m in maps]

    fig = Figure(; figure_kwargs...)
    axs = [fig[1, i] for i in 1:length(plts)]

    map(axs, plts) do ax, p
        draw!(ax, p, axis=axis)
    end

    fig
end