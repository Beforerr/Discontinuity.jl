module Discontinuity

using Arrow
using DataFrames,
    DataFramesMeta
using AlgebraOfGraphics,
    Makie
using LaTeXStrings
using Beforerr

include("naming.jl")
include("io.jl")
include("mapping.jl")
include("plot.jl")
include("formulary.jl")

export load
export plot_dist

end