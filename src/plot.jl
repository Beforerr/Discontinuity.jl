function waiting_time(time; δt=Dates.Minute(1))
    # unique and order the time
    τ = diff(time |> unique |> sort)
    return τ ./ δt
end

waiting_time(df::AbstractDataFrame, col=:time; kwargs...) = waiting_time(df[!, col]; kwargs...)

function plot_wt_pdf end
function plot_wt_pdf! end