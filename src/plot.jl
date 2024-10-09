log_tickformat = values -> [L"10^{%$(value)}" for value in values]
# axis=(xtickformat = log_tickformat,)
using AlgebraOfGraphics
using AlgebraOfGraphics: density, FigureGrid

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
    draw_kwargs = Dict()
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
    plotopts = PlotOpts(),
    kwargs...
)
    axs = [f[i, j] for j in 1:size(maps, 2) for i in 1:size(maps, 1)]
    grids = plot_dist!(axs, l, maps; kwargs...)
    fg = FigureGrid(f, grids[1])
    process_opts!(fg, axs, plotopts)
    return f
end

plot_dist!(l::Layer, maps; kwargs...) = plot_dist!(current_figure(), l, maps; kwargs...)

plot_dist(l::Layer, maps; figure=(;), kwargs...) = plot_dist!(Figure(; figure...), l, maps; kwargs...)

"""backward compatibility"""
plot_dist(l::Layer; maps=[l_log_map l_norm_log_map j_log_map j_norm_log_map], kwargs...) = plot_dist(l, maps; figure=FIGURE_KWARGS, kwargs...)

function waiting_time(time; δt=Dates.Minute(1))
    # unique and order the time
    τ = diff(time |> unique |> sort)
    return τ ./ δt
end

waiting_time(df::AbstractDataFrame, col=:time; kwargs...) = waiting_time(df[:, col])

"""Use regex to remove content within curly braces including the braces"""
format(d::Distribution) = replace(repr(d), r"\{[^}]+\}" => "")

"""
Plot the waiting time distribution of the data
"""
function plot_wt_pdf(
    τ;
    dist=Exponential,
    step=5,
    xscale=identity
)
    binedges = 0:step:maximum(τ)

    h = Hist1D(τ; binedges=binedges)
    h = normalize(h)
    d = fit(dist, τ)

    x = bincenters(h)
    y = pdf(d, x)

    f = Figure()
    ax = Axis(f[1, 1], yscale=log10, xscale=xscale, xlabel="τ (minutes)", ylabel="p(τ)")
    errorbars!(ax, h; color=:black, whiskerwidth=6)
    # plot the fit
    lines!(ax, x, y, color=:red, linewidth=2, label=format(d))

    # add the legend (fit parameters)
    axislegend(ax, loc=:upperright)

    return f
end

plot_wt_pdf(df::AbstractDataFrame; kwargs...) = plot_wt_pdf(waiting_time(df); kwargs...)