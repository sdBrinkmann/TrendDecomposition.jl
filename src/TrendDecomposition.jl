module TrendDecomposition

using LinearAlgebra
using SparseArrays
using HypothesisTests
using Statistics
using Optim

include("HPFilter.jl")
include("movingAverages.jl")
include("expSmoothing.jl")
include("l1TrendFilter.jl")

greet() = print("Hello World!")

export
    hpFilter,
    bhpFilter,
    bohlmannFilter,
    rollingAverage,
    maSeason,
    maDecompose,
    expSmoothing,
    brownLinear,
    holtLinear,
    holtWinters,
    trendL1Filter,
    trendADMM

end # module
