module TrendDecomposition

using LinearAlgebra
using SparseArrays
using HypothesisTests
using Statistics

include("HPFilter.jl")
include("movingAverages.jl")
greet() = print("Hello World!")

export HP, bHP, bohl_filter, rollingAverage, maSeason, maDecompose

end # module
