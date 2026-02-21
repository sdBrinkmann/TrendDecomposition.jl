

"""
    periodogram(Y :: Vector)

Estimates the spectral density function of a process Y with discrete spectra.

As estimation procedure the discrete fourier transform of the sample autocovariance is used.

Returns a tulpe (ω, spectra), where ω is the vector of frequency the periodogram ordinates,
spectra, are estimated for.
"""
function periodogram(Y :: Vector)
    
    T = length(Y)
    M = Int(ceil((T-1) / 2))
    y_ave = mean(Y)
    γ = zeros(Float64, T - 1);
    γ₀ = 1

    Y_cent = Y .- y_ave
    γ₀ = Y_cent' * Y_cent / T

    
    Threads.@threads for j = 1:(T-1)
        γ[j] = Y_cent[1:(T-j)]' * Y_cent[j+1:T] / T
    end
    
        
    ω = [2*i*π/T for i in 1:M]
    perio = zeros(Float64, M)

    Threads.@threads for i in 1:Int(M)
        perio[i] = 1 /
            (2*π) * (γ₀ + 2 * sum([γ[j] * cos(ω[i]*j) for j in 1:(T-1)]))
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
            (abs(1 + Φ' * [exp(-im * ω[i] * j) for j in 1:p]).^2 * 2 * π)
    end
    
   
    return (ω, spect)   
end



"""
    arSpectrum(y :: Vector; p::Int = 3, method=:burg)

Computes the spectral density of a p-th order autoregressive process; first 
the method given (:burg, :ols) is used to estimate the parameters for the time
series y and in the second step the spectral density is computed.
"""
function arSpectrum(y :: Vector; p::Int = 3, method=:burg)

    if method == :burg
        Φ, σ = arBurg(y, p)
    elseif method == :ols
        Φ, σ = arOLS(y, p)
    else
        throw(ArgumentError("Invalid method. Valid options are :burg, :ols"))
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
