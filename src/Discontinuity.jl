module Discontinuity

using DataFrames,
    DataFramesMeta
using DimensionalData
using DimensionalData: dims, TimeDim
using DimensionalData.Dimensions: Dimension
using LaTeXStrings
using LinearAlgebra
using Dates
using Distributions, FHist
using Unitful

Unitful.preferunits(u"km")

include("utils.jl")
include("naming.jl")
include("io.jl")
include("mapping.jl")
include("formulary.jl")
include("detection/detection.jl")
include("processing.jl")
include("products.jl")

export DataSet
export load, process!
export plot_dist!, plot_dist, plot_wt_pdf
export plot_B_mva!, plot_fit!, plot_fit
export waiting_time
export compute_params!, compute_Alfvenicity_params!, compute_anisotropy_params!

function plot_dist! end
function plot_dist end
function plot_fit! end
function plot_fit end
function plot_B_mva! end
function plot_wt_pdf end
function plot_wt_pdf! end

end
