module TrendDecomposition

using LinearAlgebra
using SparseArrays
using HypothesisTests
using Statistics

include("HPFilter.jl")
include("movingAverages.jl")
include("expSmoothing.jl")
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
    holtWinters
    

end # module
