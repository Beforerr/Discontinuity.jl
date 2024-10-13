module Discontinuity

using DataFrames,
    DataFramesMeta
using AlgebraOfGraphics,
    Makie
using LaTeXStrings
using LinearAlgebra
using Beforerr
using Dates
using Distributions, FHist

include("utils.jl")
include("naming.jl")
include("io.jl")
include("mapping.jl")
include("plot.jl")
include("formulary.jl")

export DataSet
export load, process!
export plot_dist, plot_wt_pdf
export waiting_time

end