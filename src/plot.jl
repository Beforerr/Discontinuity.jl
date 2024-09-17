log_tickformat = values -> [L"10^{%$(value)}" for value in values]
# axis=(xtickformat = log_tickformat,)
using AlgebraOfGraphics
using AlgebraOfGraphics: density

"""
Plot the density distribution of the thickness and current density (default) of the data layer.
"""
function plot_dist(
    data_layer;
    maps=[l_log_map l_norm_log_map j_log_map j_norm_log_map],
    axis=(yscale=log10,),
    datalimits=extrema,
    fig=missing,
    figure_kwargs=(size=(1200, 300),),
    visual=visual(Lines)
)

    plt = data_layer * density(datalimits=datalimits)
    plt = plt * visual
    plts = [plt * mapping(m) for m in maps]

    fig = ismissing(fig) ? Figure(; figure_kwargs...) : fig
    axs = [fig[i, j] for j in 1:size(maps, 2) for i in 1:size(maps, 1)]
    grids = map(axs, plts) do ax, p
        draw!(ax, p, axis=axis)
    end
    add_labels!(axs)
    pretty_legend!(fig, grids[1])

    fig
end

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