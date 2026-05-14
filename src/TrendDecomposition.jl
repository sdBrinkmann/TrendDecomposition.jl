module TrendDecomposition

using LinearAlgebra
using SparseArrays
using HypothesisTests
using Statistics
using Optim
using FFTW

include("HPFilter.jl")
include("movingAverages.jl")
include("expSmoothing.jl")
include("l1TrendFilter.jl")
include("tautString.jl")
include("arEstimation.jl")
include("spectralAnalysis.jl")
include("bandPass.jl")
include("stochastic.jl")

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
    trendADMM,
    tautADMM,
    fusedADMM,
    tautStringFit,
    
    greatestConvexMinorant,
    leastConcaveMajorant,

    arSpectrum,
    periodogram,
    arBurg,
    arOLS,
    arYuleWalker,
    arDurbinLevinson,
    autoCovariance,

    baxterKing,

    beveridgeNelson

end # module
