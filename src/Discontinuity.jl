module Discontinuity

using Arrow
using DataFrames,
    DataFramesMeta
using AlgebraOfGraphics,
    CairoMakie
using Beforerr
using LaTeXStrings


include("io.jl")
include("mapping.jl")
include("plot.jl")

export load
export plot_dist

end