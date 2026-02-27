


"""
    arBurg(y :: Vector, p :: Int; intercept::Bool = false) 

Fits an autoregressive model of order p to time series y
using Burg's method.

Returns the tulpe (Φ, σ²), where Φ is the vector of estimated coefficients
and σ² is the variance of the error terms.  
"""
function arBurg(y :: Vector, p :: Int;
                intercept::Bool = false, pacf::Bool=false) 

    n = length(y)  
    if p >= n
        throw(DomainError(p, "Order has to be lower than length of time series y"))
    elseif p <= 0
        throw(DomainError(p, "Order has to be a positve integer"))
    end

    μ = mean(y)
    X = y .- μ
    
    Θ = zeros(p, p)  # Array{Float64, 2}(undef, p, p)

    Θ[1, 1] = sum(X[2:end] .* X[1:(end-1)]) /
        sqrt(sum(X[2:end].^2) * sum(X[1:(end-1)].^2))

    W = Vector{Float64}(undef, n)
    Z = Vector{Float64}(undef, n)

    for i in 2:p
        for t in (i+1):n
            W[t] = X[t] - Θ[i-1, 1:i-1]' * X[t-1:-1:t-i+1] #X[t-i+1:t-1]
        end

          for t in 1:(n-i)
            Z[t] = X[t] - Θ[i-1, 1:i-1]' * X[t+1:t+i-1]
        end

        Θ[i, i] = 2 * sum(W[i+1:end] .* Z[1:(end-i)])  /
            (sum(W[(i+1):end].^2 .+ Z[1:(end-i)].^2))
        
        for j in 1:(i-1)
            Θ[i, j] = Θ[i-1, j] - Θ[i, i] * Θ[i-1, i-j]
        end
    end
   
    θ₀ = μ * (1 - sum(Θ[p, :]))
    Φ = Θ[p, p:-1:1]
    #σ2 = sum([y[i] - Φ' * y[i-p:i-1] - θ₀ for i in (p+1):n].^2) / (n - 2*p - 1)
    σ² = sum([X[i] - Φ' * X[i-p:i-1] for i in (p+1):n].^2) / (n - 2*p -1)
    
    if pacf == true
        if intercept == true
            return (vcat(θ₀, Θ[p, 1:p]), σ², diag(Θ))
        else        
            return (Θ[p, 1:p], σ², diag(Θ))
        end
    else
        if intercept == true
            return (vcat(θ₀, Θ[p, 1:p]), σ²)
        else        
            return (Θ[p, 1:p], σ²)
        end
    end
end



"""
    arOLS(y :: Vector, p :: Int; intercept::Bool = false) 

Fits an autoregressive model of order p to time series y
using ordinary least square (OLS).

Returns the tulpe (Φ, σ²), where Φ is the vector of estimated coefficients
and σ² is the variance of the error terms.  
"""
function arOLS(y :: Vector, p :: Int; intercept::Bool = false)

    n = length(y)

    if p >= (1/2 * n)
        throw(DomainError(p, "Order has to be lower than half of length of time series y"))
    elseif p <= 0
        throw(DomainError(p, "Order has to be a positve integer"))
    end
   
    if intercept == true
        X = hcat(ones(n-p), [y[(p+1-i):(end-i)] for i in 1:p]...)
        Φ = inv(X'X) * X' * y[p+1:end]
        σ2 = sum([y[i] - Φ' * vcat(1, y[i-1:-1:i-p]) for i in (p+1):n].^2) / (n - 2*p -1)
    else
        μ = mean(y)
        Y = y .- μ
        X = hcat([Y[(p+1-i):(end-i)] for i in 1:p]...)
        Φ = inv(X'X) * X' * Y[p+1:end]
        #Φ = X'X \ X' * Y[p+1:end]
        #param = reverse(Φ')
        σ² = sum([Y[i] - Φ' * Y[i-1:-1:i-p] for i in (p+1):n].^2) / (n - 2*p -1)
    end

    return (Φ, σ²)
    
end


"""
    arYuleWalker(y :: Vector, p :: Int)

Fits an autoregressive model of order p to time series y
using the Yule-Walker estimators.

Returns the tulpe (Φ, σ²), where Φ is the vector of estimated coefficients
and σ² is the variance of the error terms.  
"""
function arYuleWalker(y :: Vector, p :: Int)

    T = length(y)

    if p >= (1/2 * T)
        throw(DomainError(p, "Order has to be lower than half of length of time series y"))
    elseif p <= 0
        throw(DomainError(p, "Order has to be a positve integer"))
    end
    μ = mean(y)
    X = y .- μ

    γ₀ = X' * X / T
    γ = autoCovariance(X, T, maxCov=p);

    function indexCov(index :: Int)
        if index == 0
            return γ₀
        elseif index > 0
            return γ[index]
        else
            return γ[abs(index)]
        end
    end
  
    Γ = zeros(p, p)

    for j in 1:p
        for i in 1:p
            Γ[i, j] = indexCov(i - j)
        end
    end
  
    Φ = inv(Γ) * γ[1:p]
    σ² =  γ₀ - Φ' * γ[1:p]    
    #σ2 = sum([Y[i] - Φ' * Y[i-1:-1:i-p] for i in (p+1):T].^2) / (T - 2*p -1)

    return (Φ, σ²)
end


"""
    arYuleWalker(y :: Vector, p :: Int)

Fits an autoregressive model of order p to time series y
using the Durbin-Levinson recursion.

Returns the tulpe (Φ, σ²), where Φ is the vector of estimated coefficients
and σ² is the variance of the error terms.  
"""
function arDurbinLevinson(y :: Vector, p :: Int;
                          intercept::Bool = false, pacf::Bool = false) 

    n = length(y)
    
    if p >= n
        throw(DomainError(p, "Order has to be lower than length of time series y"))
    elseif p <= 0
        throw(DomainError(p, "Order has to be a positve integer"))
    end
 
    μ = mean(y)
    X = y .- μ

    γ₀ = X' * X / n
    γ = autoCovariance(X, n, maxCov=p)
    
    Θ = zeros(p, p)  # Array{Float64, 2}(undef, p, p)

    Θ[1, 1] = γ[1] / γ₀
    σ = γ₀ * (1 - Θ[1, 1].^2)
    #yₙ * (1 - Θ[1,1])
    
    for i in 2:p
        Θ[i, i] = (γ[i] - sum([Θ[i-1, j] * γ[i-j] for j in 1:(i-1)])) / σ      
        for j in 1:(i-1)
            Θ[i, j] = Θ[i-1, j] - Θ[i, i] * Θ[i-1, i-j]
        end
        σ = σ * (1 - Θ[i, i].^2)        
    end
   
    θ₀ = μ * (1 - sum(Θ[p, :]))
    Φ = Θ[p, p:-1:1]
    #σ2 = sum([y[i] - Φ' * y[i-p:i-1] - θ₀ for i in (p+1):n].^2) / (n - 2*p - 1)
    #σ² = sum([X[i] - Φ' * X[i-p:i-1] for i in (p+1):n].^2) / (n - 2*p -1)

    if pacf == true
        if intercept == true
            return (vcat(θ₀, Θ[p, 1:p]), σ, diag(Θ))
        else        
            return (Θ[p, 1:p], σ, diag(Θ))
        end
    else
        if intercept == true
            return (vcat(θ₀, Θ[p, 1:p]), σ)
        else        
            return (Θ[p, 1:p], σ)
        end
    end
end
