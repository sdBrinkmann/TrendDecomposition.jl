# ARMA.jl 


spL(T, q) = spdiagm(T, T, -q => fill(1., T-q))


function colInv!(a :: Vector, α :: Vector, q :: Int, T :: Int)
    a[1] = 1.0
    for t in 2:q
        #acc = 0.0
        a[t] = 0.
        for s in 1:t-1
            @inbounds a[t] -= α[s] * a[t-s, 1]
        end
    end
    for t in (q+1):T
        a[t] = 0.
        for s in 1:q
            @inbounds a[t] -= α[s] * a[t-s, 1]
        end
    end
end




function mulB!(v::Vector, β::Vector, Ay::Vector, p :: Int, T::Int)
    a = [1.0; β]

    for idx in 1:p
        v[idx] = a[1:idx]' * Ay[idx:-1:1]
    end
    
    Threads.@threads for idx in p+1:T
        v[idx] = 0.0 
        for j in 1:p+1
            @inbounds v[idx] += a[j] * Ay[idx-j+1]
        end
    end
end

    
function mulA!(v::Vector, a::Vector, y::Vector, T::Int)
    Threads.@threads for idx in 1:T
        v[idx] = 0.0
        for j in 1:idx
            @inbounds v[idx] += a[j] * y[idx-j+1]
        end
    end
end

function mulL(a::Vector, b::Vector, l :: Int, T::Int)
    acc = 0.0
    for i in 1:(T-l)
        acc += a[i+l] * b[i]
    end    
    return acc
end

function mulLL(a::Vector, b::Vector, l::Int, k::Int, T::Int)
    acc = 0.0
    if l >= k
        for i in 1:(T-l)
            @inbounds acc += a[i] * b[i+l-k]
        end
    else
        dif = k - l
        for i in 1:(T-k)
            @inbounds acc += a[i+dif] * b[i]
        end
    end
    return acc   
end

"""
    maNR(y :: Vector, q::Int = 1; iter::Int = 5, α = fill(0.2, q)) 

Estimates a MA(q) model using the maximum likelihood function and
the Newton-Raphson procedure.

Returns a tuple, (α, σ) where α is the (qx1) vector containing the MA  coeffients,
and σ is the estimated variance. 
"""
function maNR(y :: Vector, q::Int = 1; iter::Int = 5, α = fill(0.2, q))   
    T :: Int = length(y)
    #α :: Vector{Float64} = fill(0.4, q)
    q_vec = Vector{Float64}(undef, q)

    Φ = zeros(Float64, q, q)
    a = Vector{Float64}(undef, T) 
    Ay = Vector{Float64}(undef, T)
 
    V1 = Vector{Float64}(undef, T)

    for i in 1:iter 

        colInv!(a, α, q, T) 
        mulA!(Ay, a, y, T)
        mulA!(V1, a, Ay, T)
        
        Threads.@threads for j in 1:q
            q_vec[j] =  mulL(Ay, V1, j, T) 
            for g in 1:q 
                Φ[g,j] = mulLL(V1, V1, g, j, T) 
            end
        end
        
        sol = Φ \ q_vec
        α = α + sol
    end

    colInv!(a, α, q, T) 
    mulA!(Ay, a, y, T) 

    σ :: Float64 = Ay' * Ay / T
    
    return (α, σ)
end


"""
    armaNR(y :: Vector, p::Int = 1, q::Int = 1; iter::Int = 5, α = fill(0.4, q), β=[-1.])

Estimates an ARMA(p, q) model using the maximum likelihood function and
the Newton-Raphson procedure.

Returns a tuple, (β, α, σ) where β and α are the vectors containing the p AR and q MA coefficients respectively, and σ is the estimated variance. 
"""
function armaNR(y :: Vector, p::Int = 1, q::Int = 1; iter::Int = 5, α = fill(0.4, q), β=[-1.])
    
    T :: Int = length(y)
    #α :: Vector{Float64} = fill(0.0, q)
    #β :: Vector{Float64} = fill(0.0, p)

    if β == [-1.]
        β::Vector{Float64}, _ = arYuleWalker(y, p)
        β = -β
    end

    qvec = Vector{Float64}(undef, q)
    pvec = Vector{Float64}(undef, p)

    Φ = zeros(Float64, q, q)
    Ψ = zeros(Float64, p, p)
    Ω = zeros(Float64, q, p)
    
    a = Vector{Float64}(undef, T) 
    Ay = Vector{Float64}(undef, T) 
    v = Vector{Float64}(undef, T)
    V1 = Vector{Float64}(undef, T)

    for i in 1:iter
       
        colInv!(a, α, q, T)
        mulA!(Ay, a, y, T)
        mulB!(v, β, Ay, p, T)
        mulA!(V1, a, v, T)
        
        for g in 1:q
            qvec[g] = mulL(v, V1, g, T) 
        end

        for k in 1:p
            pvec[k] = - mulL(v, Ay, k, T)
        end
     
         Threads.@threads for f in 1:q
             for g in 1:q
                 Φ[g, f] = mulLL(V1, V1, g, f, T)
             end
        end

        Threads.@threads for l in 1:p
            for g in 1:q
                Ω[g, l] = - mulLL(V1, Ay, g, l, T) 
            end
        end
        
        mulA!(V1, a, y, T)
        
        Threads.@threads for k in 1:p
            for l in 1:p
                Ψ[k, l] = mulLL(V1, Ay, k, l, T) 
            end
        end 
            
        sol = [Φ Ω; Ω' Ψ] \ [qvec; pvec]
        
        α = α .+ sol[1:q]
        β = β .+ sol[q+1:end]      
    end
    
    colInv!(a, α, q, T)
    mulA!(Ay, a, y, T)
    mulB!(v, β, Ay, p, T)
    σ = v' * v / T
    return (-β, α, σ)
end



function likelihoodARMA(y :: Vector, ϕ :: Vector, θ :: Vector, σ :: Float64)

    T :: Int = length(y)
    a = Vector{Float64}(undef, T) 
    Ay = Vector{Float64}(undef, T)
    v = Vector{Float64}(undef, T)

    q = length(θ)
    p = length(ϕ)
    
    colInv!(a, θ, q, T)
    mulA!(Ay, a, y, T)
    mulB!(v, ϕ, Ay, p, T)

    L = - 0.5 * T * log(2*π)  - 0.5 * T * log(σ)  -  v' * v / (σ * 2)
  
    return L
end


function likelihoodMA(y :: Vector, θ :: Vector, σ :: Float64)

    T :: Int = length(y)
    a = Vector{Float64}(undef, T) 
    Ay = Vector{Float64}(undef, T)
    q = length(θ)

    
    colInv!(a, θ, q, T)
    mulA!(Ay, a, y, T)

    L = - 0.5 * T * log(2*π)  - 0.5 * T * log(σ)  -  Ay' * Ay / (σ * 2)
  
    return L
end
