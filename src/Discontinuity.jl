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
using Unitful

Unitful.preferunits(u"km")

include("utils.jl")
include("naming.jl")
include("io.jl")
include("mapping.jl")
include("plot.jl")
include("formulary.jl")
include("detection/detection.jl")
include("processing.jl")
include("products.jl")

export DataSet
export load, process!
export plot_dist, plot_wt_pdf
export waiting_time

end
