
"""
    beveridgeNelson(X :: Vector, p :: Int; method::Symbol=:Newbold, estimator::Symbol=:OLS)

Computes the Beverdige-Nelson decomposition of a non-stationary time series X using a
autoregressive model for estimation with p-lags, as estimators :OLS, :burg, :yuleWalker
can be used.

Either uses Newbold's (1990) or Miller's (1987) method of computation. 

Returns a (Nx1) vector of the trend component.
"""
function beveridgeNelson(X :: Vector, p :: Int; method::Symbol=:Newbold, estimator::Symbol=:OLS)

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
    
    Ω = 1 / (1 - sum(ϕ[2:end]))
    m = sum([j * ϕ[j+1] * ϕ[1] for j in 1:p])
    w = vcat(1, -ϕ[2:end])
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
    dX_m = dX .- ϕ[1]
    c = [W' * dX_m[i-1:-1:i-p] for i in (p+2):(n)]
    return vcat(repeat([NaN], outer=p+1), c + X[p+2:end])
end


