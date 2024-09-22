log_tickformat = values -> [L"10^{%$(value)}" for value in values]
# axis=(xtickformat = log_tickformat,)
using AlgebraOfGraphics
using AlgebraOfGraphics: density

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
    visual=visual(Lines)
)
    plt = l * density(datalimits=datalimits)
    plt = plt * visual
    plts = [plt * mapping(m) for m in maps]
    return map(axs, plts) do ax, p
        draw!(ax, p, axis=axis)
    end
end


function plot_dist!(
    f::Figure,
    l::Layer,
    maps;
    add_labels = true,
    kwargs...
)
    axs = [f[i, j] for j in 1:size(maps, 2) for i in 1:size(maps, 1)]
    grids = plot_dist!(axs, l, maps; kwargs...)
    add_labels || add_labels!(axs)
    pretty_legend!(f, grids[1])
    return f
end

"""backward compatibility"""
plot_dist!(f::Figure,l::Layer; maps=[l_log_map l_norm_log_map j_log_map j_norm_log_map], kwargs...) = plot_dist!(f, l, maps; kwargs...)

plot_dist!(l::Layer; kwargs...) = plot_dist!(current_figure(), l; kwargs...)
plot_dist(l::Layer; figure=FIGURE_KWARGS, kwargs...) = plot_dist!(Figure(; figure...), l; kwargs...)

function waiting_time(time; δt=Dates.Minute(1))
    # unique and order the time
    τ = diff(time |> unique |> sort)
    return τ ./ δt
end

waiting_time(df::AbstractDataFrame; col=:time, kwargs...) = waiting_time(df[:, col])

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