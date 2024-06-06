log_tickformat = values -> [L"10^{%$(value)}" for value in values]
# axis=(xtickformat = log_tickformat,)
using AlgebraOfGraphics
using AlgebraOfGraphics: density

function hideylabels(la::Axis)
    la.ylabelvisible = false
end

function hideylabels!(fgs)
    if length(fgs) > 1
        [ hideylabels.(fg) for fg in fgs[2:end] ]
    end
end

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
    fgs = [draw!(ax, p, axis=axis) for (ax, p) in zip(axs, plts)]

    # hide extra y labels
    hideylabels!(fgs)
    add_labels!(axs)
    pretty_legend!(fig, fgs[1])
    fig
end

for sym in [:hideylabels,]
    @eval function $sym(ae::AxisEntries; kwargs...)
        axis = ae.axis
        AlgebraOfGraphics.isaxis2d(axis) && $sym(axis; kwargs...)
    end
end