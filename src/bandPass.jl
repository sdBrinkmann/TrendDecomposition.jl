

"""
    baxterKing(x :: Vector, pl::Real, pu::Real, K::Int)

Computes the Baxter-King band pass filter; pu and pl determine the
lower and upper cutoff frequency, respectively, in terms of the periodicity.
K denotes the number of included lags.

Returns a (N-2K x 1) vector of the filtered cycle components. 
"""
function baxterKing(x :: Vector, pl::Real, pu::Real, K::Int)

    if pl < 2
        throw(DomainError(pl, "Lower bound has to be greater equal than 2"))
    elseif pu < pl
        throw(DomainError(pu, "Upper bound has to be greater than lower bound"))
    end
    
    a1 = 2*π / pu
    a2 = 2*π / pl
    
    b = zeros(K)
    b0 = (a2 - a1) / π
    for i in 1:K
        b[i] = (sin(i * a2) - sin(i * a1)) / (π * i)
    end
    W = vcat(reverse(b), b0, b)
    term = sum(W) / (2 * K + 1)
    W = W .- term
    return rollingAverage(x, 2*K + 1, W, discard=true)
end
