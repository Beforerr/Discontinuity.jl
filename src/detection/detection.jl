using DataFrames, DataFramesMeta
using DimensionalData
using DimensionalData: dims
using DimensionalData.Dimensions: Dimension
using LinearAlgebra
using Statistics
using NaNStatistics
using StaticArrays
using OhMyThreads

export detect_variance, compute_indices, filter_indices

include("utils.jl")
include("variance.jl")