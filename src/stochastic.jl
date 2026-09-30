
"""
    beveridgeNelson(X :: Vector, p :: Int; method::Symbol=:Newbold, estimator::Symbol=:OLS)

Computes the Beverdige-Nelson decomposition of a non-stationary time series X using an
autoregressive model for estimation with p-lags. As estimators :OLS, :burg, :yuleWalker
can be used.

Either uses Newbold's (1990) or Miller's (1987) method of computation. 

Returns a (Nx1) vector of the trend component.
"""
function beveridgeNelson(X :: Vector, p::Int; method::Symbol=:Newbold, estimator::Symbol=:OLS)

    if method == :Newbold
        return bnNewbold(X, p, estimator=estimator)
    elseif method == :Miller
        return bnMiller(X, p, estimator=estimator)
    else
        throw(ArugmentError("Invalid Estimator. Valid options are :Newbold or :Miller"))
    end
end


function bnMiller(X :: Vector, p :: Int; estimator::Symbol=:OLS)
    n = length(X)
    dX = X[2:end] .- X[1:(n-1)]
    
    if estimator == :burg
        ϕ, _ = arBurg(dX, p)
    elseif estimator == :OLS
        ϕ, _ = arOLS(dX, p)
    elseif estimator == :yuleWalker
        ϕ, _ = arYuleWalker(dX, p)
    elseif estimator == :durbinLevinson
        ϕ, _ = arDurbinLevinson(dX, p)
    else
        throw(ArgumentError("Invalid estimator. Valid options are :burg, :OLS, :yuleWalker"))
    end
    
    Ω = 1 / (1 - sum(ϕ))
    w = vcat(1, -ϕ) .* Ω
    τ = [w' * X[i:-1:i-p] for i in (p+1):n]
    return vcat(repeat([NaN], outer=p), τ)
end


function bnNewbold(X :: Vector, p :: Int; estimator::Symbol=:OLS)
    n = length(X)
    dX = X[2:end] .- X[1:(n-1)]
       
    if estimator == :burg
        ϕ, _ = arBurg(dX, p)
    elseif estimator == :OLS
        ϕ, _ = arOLS(dX, p)
    elseif estimator == :yuleWalker
        ϕ, _ = arYuleWalker(dX, p)
    elseif estimator == :durbinLevinson
        ϕ, _ = arDurbinLevinson(dX, p)
    else
        throw(ArgumentError("Invalid estimator. Valid options are :burg, :OLS, :yuleWalker"))
    end
    
    Ω = 1 / (1 - sum(ϕ))
    μ = mean(dX)
    m = sum([j * ϕ[j] * μ for j in 1:p])
    w = vcat(1, -ϕ)
    τ =  [w' * X[i:-1:i-p] - m for i in (p+1):n]
    return vcat(repeat([NaN], outer=p), τ .* Ω)
end


function bnNewbold2(X :: Vector, p :: Int; estimator::Symbol=:OLS)
    n = length(X)
    dX = X[2:end] .- X[1:(n-1)]

    
    if estimator == :burg
        ϕ, _ = arBurg(dX, p, intercept=true)
    elseif estimator == :OLS
        ϕ, _ = arOLS(dX, p, intercept=true)
    elseif estimator == :yuleWalker
        ϕ, _ = arYuleWalker(dX, p, intercept=true)
    elseif estimator == :durbinLevinson
        ϕ, _ = arDurbinLevinson(dX, p, intercept=true)
    else
        throw(ArgumentError("Invalid estimator. Valid options are :burg, :OLS, :yuleWalker"))
    end

    e = vcat(1, zeros(p-1))
    A = [reshape(ϕ[2:end], 1, p); Diagonal(ones(p))]
    A = A[1:p, :]
    W = vec(e' * inv(I - A) * A)
    μ = mean(dX)
    dX_m = dX .- μ 
    c = [W' * dX_m[i-1:-1:i-p] for i in (p+1):(n)]
    return vcat(repeat([NaN], outer=p), c + X[p+1:end])
end

"""
    armaForecast(y :: Vector, ϕ :: Vector, θ :: Vector; h::Int = 1)

Forecast time series y (with length T) as an ARMA(p, q) process with parameter vectors
ϕ and θ. The forecast horizon is given by h.

Returns (hx1) vector of forecast values for time points T+1, ..., T+h.
"""
function armaForecast(y :: Vector, ϕ :: Vector, θ :: Vector; h::Int = 1)
    if h < 1
        throw(ArgumentError("forecast horizion h has to be greater than 0"))
    end
    p = length(ϕ)
    q = length(θ)
    T = length(y)
    ϵ = zeros(Float64, T)
    Y = [y; zeros(Float64, h)]

    ict = mean(y) * (1 - sum(ϕ))

    if p >= q
        for i in (p+1):T
            ϵ[i] = y[i] - ict - ϕ' * y[i-1:-1:i-p] - θ' * ϵ[i-1:-1:i-q] 
        end
    else
        for i in (p+1):q
            ϵ[i] = y[i] - ict - ϕ' * y[i-1:-1:i-p] - θ[1:i-1]' * ϵ[i-1:-1:1] 
        end 
        
        for i in (q+1):T
            ϵ[i] = y[i] - ict - ϕ' * y[i-1:-1:i-p] - θ' * ϵ[i-1:-1:i-q] 
        end
    end

    for (off ,i) in enumerate(T+1:T+min(q, h))
        Y[i] = ict + ϕ' * Y[i-1:-1:i-p] + θ[1+off-1:q]' * ϵ[i-off:-1:i-q]        
    end

    if h > q
        for i in T+1+q:T+h
            Y[i] = ϕ' * Y[i-1:-1:i-p]  + ict  
        end
    end
    
    return Y[T+1:T+h]
end


"""
    maForecast(y :: Vector, θ :: Vector; h::Int = 1)

Forecast time series y (with length T) as an MA(q) process with parameter vector θ.
The forecast horizon is given by h.

Returns (max(h, q)x1) vector of forecast values for time points T+1, ..., T+max(h, q).
"""
function maForecast(y :: Vector, θ :: Vector; h::Int = 1)
    if h < 1
        throw(ArgumentError("forecast horizion h has to be greater than 0"))
    end

    q = length(θ)
    T = length(y)
    ϵ = zeros(Float64, T)

    Y = zeros(Float64, min(h, q))

    μ = mean(y)
    ϵ[1] = y[1] - μ
    
    for i in 2:q
        ϵ[i] = y[i] - μ - θ[1:i-1]' * ϵ[i-1:-1:1] 
    end

    for i in (q+1):T
        ϵ[i] = y[i] - μ - θ' * ϵ[i-1:-1:i-q] 
    end

    for (off ,i) in enumerate(T+1:T+min(q, h))
        Y[off] = μ + θ[1+off-1:q]' * ϵ[i-off:-1:i-q]        
    end
  
    return Y
end

function bnNewbold(X :: Vector, p :: Int, q :: Int)
    
    n = length(X)
    dX = X[2:end] .- X[1:(n-1)]
       
    ϕ, θ , _ = armaNR(dX, p, q)
    μ = mean(dX)
     

    e = vcat(1, zeros(p-1))
    A = [reshape(ϕ, 1, p); Diagonal(ones(p))]
    A = A[1:p, :]
    W = vec(e' * inv(I - A) * A)
    μ = mean(dX)
    dX_m = dX .- μ 
    c = zeros(Float64, n - p)
    fore = zeros(Float64, q)
    
    for i in (p+1):n
        fore = armaForecast(dX[1:i-1], ϕ, θ, h = q) .- μ
        if p <= q
            c[i-p] = W' * fore[q:-1:q-p+1] + sum(fore)
        else
            c[i-p] = W[1:q]' * fore[q:-1:1] + W[q+1:end]' * dX_m[i-1:-1:i-p+q]  + sum(fore)
        end
    end
    
    return vcat(repeat([NaN], outer=p), c + X[p+1:end])
end


function bnNewbold(X :: Vector, ϕ :: Vector, θ :: Vector)
    
    n = length(X)
    dX = X[2:end] .- X[1:(n-1)]

    p = length(ϕ)
    q = length(θ)

    μ = mean(dX)
     

    e = vcat(1, zeros(p-1))
    A = [reshape(ϕ, 1, p); Diagonal(ones(p))]
    A = A[1:p, :]
    W = vec(e' * inv(I - A) * A)

    μ = mean(dX)
    dX_m = dX .- μ 
    c = zeros(Float64, n - p)
    fore = zeros(Float64, q)
    
    for i in (p+1):n
        fore = armaForecast(dX[1:i-1], ϕ, θ, h = q) .- μ
        if p <= q
            c[i-p] = W' * fore[q:-1:q-p+1] + sum(fore)
        else
            c[i-p] = W[1:q]' * fore[q:-1:1] + W[q+1:end]' * dX_m[i-1:-1:i-p+q]  + sum(fore)
        end
    end
    
    return vcat(repeat([NaN], outer=p), c + X[p+1:end])
end



"""
    beveridgeNelson(X :: Vector, p :: Int, q :: Int)

Computes the Beverdige-Nelson decomposition of a non-stationary time series X using an
ARMA(p, q) model. 

Uses Newbold's (1990) method of computation. 

Returns a (Nx1) vector of the trend component.
"""
function beveridgeNelson(X :: Vector, p :: Int, q :: Int)

    if q > 0 && p > 0
        return bnNewbold(X, p, q)
    elseif p > 0
        return beveridgeNelson(X, p)
    elseif q > 0
        n = length(X)
        dX = X[2:end] .- X[1:(n-1)]
        θ , _ = maNR(dX, q)
        μ = mean(dX)
        
        ma_term = [sum(maForecast(dX[1:i], θ, h = q) .- μ) for i in q:(n-1)]
        return vcat(repeat([NaN], outer=q), ma_term .+ X[q+1:end] )
    else
        throw(ArugmentError("q and p have to be non-negative Integers"))
    end
end



"""
    beveridgeNelson(X :: Vector, ϕ :: Vector, θ :: Vector)

Computes the Beverdige-Nelson decomposition of a non-stationary time series X using an
ARMA(p, q) model of ΔX given the AR parameter ϕ and MA parameter θ. 

Uses Newbold's (1990) method of computation. 

Returns a (Nx1) vector of the trend component.
"""
function beveridgeNelson(X :: Vector, ϕ :: Vector, θ :: Vector)

    return bnNewbold(X, ϕ, θ)
end
