module Discontinuity

using Arrow
using DataFrames,
    DataFramesMeta
using AlgebraOfGraphics,
    Makie
using LaTeXStrings
using Beforerr
using Dates
using Distributions, FHist

include("naming.jl")
include("io.jl")
include("mapping.jl")
include("plot.jl")
include("formulary.jl")

export DataSet
export load
export plot_dist, plot_wt_pdf
export waiting_time

end