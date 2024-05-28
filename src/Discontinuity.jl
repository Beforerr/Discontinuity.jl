module Discontinuity

using AlgebraOfGraphics,
    CairoMakie
using AlgebraOfGraphics: density
using beforerr
using LaTeXStrings

include("mapping.jl")
include("plot.jl")
include("io.jl")

export load

end