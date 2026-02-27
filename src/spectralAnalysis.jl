

"""
    autoCovariance(X :: Vector, T::Int; maxCov::Int=T-1)
    
Calculates the autocovariance of time series X
"""
function autoCovariance(X :: Vector, T::Int; maxCov::Int=T-1)
    
    #T = length(X)   
    γ = Array{Float64, 1}(undef, maxCov)
    Threads.@threads for j = 1:maxCov
        acc = 0
        for i in 1:(T-j)
            @inbounds acc += X[i] * X[i+j]
            @inbounds γ[j] = acc / T
        end
    end
    return γ
end


"""
    periodogram(Y :: Vector; trunc::Int=-1)

Estimates the spectral density function of a process Y with discrete spectra.

As estimation procedure the discrete fourier transform of the sample autocovariance is used. A truncation point can be given with trunc, so that a truncated periodogram
is estimated. 

Returns a tulpe (ω, spectra), where ω is the vector of frequency the periodogram ordinates,
spectra, are estimated for.
"""
function periodogram(Y :: Vector; trunc::Int=-1)
    
    T = length(Y)
    M = Int(ceil((T-1) / 2))
    μ = mean(Y)
    
    if trunc >= T
        throw(DomainError(trunc, "trunc has to be smaller than length of vector Y"))
    elseif trunc == -1
        trunc = T - 1
    end
    
    X = Y .- μ
    γ₀ = X' * X / T
    γ = autoCovariance(X, T, maxCov=trunc)
    
    ω = [2*i*π/T for i in 1:M]
    perio = zeros(Float64, M)

    Threads.@threads for i in 1:Int(M)
        acc = 0
        for j in 1:trunc
            @inbounds acc += γ[j] * cos(ω[i] * j)
        end
        perio[i] = 1 / (2*π) * (γ₀ + 2 * acc )
    end

    return (ω, perio)
end



"""
    arSpectrum(Φ :: Vector, σ::Float64; T = 100)

Computes the spectral density of a p-th order autoregressive process, given the
parameter vector Φ and the variance σ. 

"""
function arSpectrum(Φ :: Vector, σ::Float64; T = 100)

    #T = length(y)
    M = Int(ceil((T-1) / 2))
    ω = [2*i*π/T for i in 1:M]

    #Φ = -1 * Φ
    p = length(Φ)
    
    spect = zeros(M)
         
    for i in 1:M 
        spect[i] =  σ /
            (abs(1 - Φ' * [exp(-im * ω[i] * j) for j in 1:p]).^2 * 2 * π)
    end
      
    return (ω, spect)   
end



"""
    arSpectrum(y :: Vector; p::Int = 3, method=:burg)

Computes the spectral density of a p-th order autoregressive process; first 
the method given (:burg, :ols, :yuleWalker, :durbinLevinson) is used to estimate
the parameters for the time series y and in the second step the spectral density is computed.
"""
function arSpectrum(y :: Vector; p::Int = 3, method=:burg)

    if method == :burg
        Φ, σ = arBurg(y, p)
    elseif method == :ols
        Φ, σ = arOLS(y, p)
    elseif method == :yuleWalker
        Φ, σ = arYuleWalker(y, p)
    elseif method == :durbinLevinson
        Φ, σ = arDurbinLevinson(y, p)
    else
        throw(ArgumentError("Invalid method. Valid options are :burg, :ols, :yuleWalker"))
    end
              
    T = length(y)
    M = Int(ceil((T-1) / 2))
    ω = [2*i*π/T for i in 1:M]

    spect = zeros(M)    
    
    for i in 1:M
        spect[i] =  σ /
            (abs(1 - Φ' * [exp(-im * ω[i] * j) for j in 1:p]).^2 * 2 * π)
    end

    return (ω, spect)   
end
