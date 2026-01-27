



"""
    expSmoothing(y :: Vector, λ :: Real; startValue::Bool = true, v_0::Real = 0.0)

Simple exponential smoothing with smoothing factor λ. If startValue equals true,
v_0 will be used as initial value, otherwise the first element of y (y[1])
will be selected as starting value.
Using y[1] as starting value will result in a different mathematical function.   
"""
function expSmoothing(y :: Vector, λ :: Real; startValue::Bool = true, v_0::Real = 0.0)
    @assert λ >= 0 && λ <= 1
    wma = zeros(length(y))
    if startValue == true
        wma[1] = (1-λ) * v_0 + λ * y[1]
    else
        wma[1] = y[1]
    end

    for t in 2:length(y)
        wma[t] = (1-λ) * wma[t-1] + λ * y[t]
    end
    return wma
end


"""
    brownLinear(y :: Vector, λ :: Real; h::Int = 0, startValues::Tuple=(), φ::Real = 1.0)

Double exponential smoothing with smoothing factor λ and damping factor
φ (\\varphi) given vector y (N x 1) and h forecast periods.

Start values of level and slope can be given as tuple (l₀, b₀),
otherwise output values start at period 3.

Returns a tuple with (f, D), where f is a (N+h) x 1 vector with the
forecast values and D is a (N x 2) matrix with the level and slope as columns.

"""
function brownLinear(y :: Vector, λ :: Real; h::Int = 0, startValues::Tuple=(), φ::Real = 1.0)
    T = length(y)
    F = zeros(T+h)
    lvl = zeros(T)
    b = zeros(T)
    start :: Int = 3
    
    if startValues == ()
        b[1] = lvl[1] = NaN
        b[2] = y[2] - y[1]
        lvl[2] = y[2]
        F[1:2] .= NaN
    else
        F[1] = startValues[1] + startValues[2] * φ
        lvl[1] = startValues[1] + startValues[2] * φ + (1. - λ^2) * (y[1] - F[1])
        b[1] =  startValues[2] * φ + (1. - λ)^2 * (y[1] - F[1])
        start = 2
    end

    for t in start:T
        F[t] = lvl[t-1] + b[t-1] * φ
        lvl[t] = lvl[t-1] + b[t-1] * φ + (1. - λ^2) * (y[t] - F[t])
        b[t] = b[t-1] * φ + (1. - λ)^2 * (y[t] - F[t])
    end

    if φ == 1.0
        for j in 1:h
            F[T + j] = lvl[T] + b[T] * j
        end
    else
        for j in 1:h
            F[T + j] = lvl[T] + b[T] * sum([φ^i for i in 1:j])
        end
    end

    return (F, hcat(lvl, b))
end

"""
    brownLinear(y :: Vector; h::Int = 0, startValues::Tuple=(), φ::Union{Real, Nothing} = nothing)

Double exponential smoothing with damping factor φ (\\varphi) given vector y (N x 1)
and h forecast periods. The smoothing parameter λ is estimated by minimizing the loss
function using Optim.jl. For φ equal nothing, its value is also determined by the optimization algorithm.

Start values of level and slope can be given as tuple (l₀, b₀),
otherwise output values start at period 3.

Returns a tuple with (f, D, p), where f is a (N+h) x 1 vector with the
forecast values and D is a (N x 2) matrix with the level and slope as columns,
and p is a vector containing the estimated parameter(s)

"""
function brownLinear(y :: Vector; h::Int = 0, startValues::Tuple=(), φ::Union{Real, Nothing} = nothing)
      start :: Int = (startValues == () ? 3 : 1)
    if φ == nothing
        function optPhi(λ)
            f, _ = brownLinear(y, λ[1], startValues=startValues, φ = λ[2])
            return mse(y[start:end], f[start:end])
        end
        lower_bound = zeros(2)
        upper_bound = ones(2)
        res = optimize(optPhi, lower_bound, upper_bound, [.8, .8])
        λ = Optim.minimizer(res)
        f, D = brownLinear(y, λ[1], startValues=startValues, h=h, φ = λ[2])
        return (forecast=f, model=D, parameters=λ)
    else
        function opt(λ)
            f, _ = brownLinear(y, λ, startValues=startValues, φ=φ)
            return mse(y[start:end], f[start:end])
        end
        res = optimize(opt, 0., 1.)
        λ = Optim.minimizer(res)
        f, D = brownLinear(y, λ, startValues=startValues, h=h, φ = φ)
        return (forecast=f, model=D, parameter=λ)
    end 
end



"""
    holtLinear(y :: Vector, λ₁ :: Real, λ₂ :: Real;
                    startValues::Tuple=(), h::Int = 0, φ::Real = 1.0)

Holt's linear trend method with smoothing factors λ₁ (level) and λ₂ (slope)
and damping factor φ (\\varphi) given vector y (N x 1) and h forecast periods.

Start values of level and slope can be given as tuple (l₀, b₀),
otherwise output values start at period 3.

Returns a tuple with (f, D), where f is a ((N+h) x 1) vector with the
forecast values and D is a (N x 2) matrix with the level and slope as columns.

"""
function holtLinear(y :: Vector, λ₁ :: Real, λ₂ :: Real;
                    startValues::Tuple=(), h::Int = 0, φ::Real = 1.0)
    T = length(y)
    F = zeros(T + h)
    lvl = zeros(T)
    b = zeros(T)

    start :: Int = 3
    
    if startValues == ()
        b[1] = lvl[1] = NaN
        b[2] = y[2] - y[1]
        lvl[2] = y[2]
        F[1:2] .= NaN
    else
        lvl[1] = λ₁*y[1] + (1- λ₁)*(startValues[1] + startValues[2] * φ)
        b[1] = λ₂ * (lvl[1] - startValues[1]) + (1 - λ₂) * startValues[2] * φ
        F[1] = startValues[1] + startValues[2] * φ
        start = 2
    end

    for t in start:T
        lvl[t] = λ₁*y[t] + (1- λ₁)*(lvl[t-1] + b[t-1] * φ)
        b[t] = λ₂ * (lvl[t] - lvl[t-1]) + (1 - λ₂) * b[t-1] * φ
        F[t] = lvl[t-1] + b[t-1] * φ
    end

    if φ == 1.0
        for j in 1:h
            F[T + j] = lvl[T] + b[T] * j
        end
    else
        for j in 1:h
            F[T + j] = lvl[T] + b[T] * sum([φ^i for i in 1:j])
        end
    end

    return (F, hcat(lvl, b))
end


"""
    holtLinear(y :: Vector;
                    startValues::Tuple=(), h::Int = 0, φ::Union{Real, Nothing} = nothing)

Holt's linear trend method with damping factor φ (\\varphi) given vector y (N x 1)
and h forecast periods.

The smoothing parameters are estimated by minimizing the loss
function using Optim.jl. For φ equal nothing, its value is also determined by the optimization algorithm.

Start values of level and slope can be given as tuple (l₀, b₀),
otherwise output values start at period 3.

Returns a tuple with (f, D, p), where f is a ((N+h) x 1) vector with the
forecast values and D is a (N x 2) matrix with the level and slope as columns,
and p is a vector containing the estimated parameter(s).

"""
function holtLinear(y :: Vector;
                    startValues::Tuple=(), h::Int = 0, φ::Union{Real, Nothing} = nothing)
    
    start :: Int = (startValues == () ? 3 : 1)
    if φ == nothing
        function optPhi(λ)
            f, _ = holtLinear(y, λ[1], λ[2], startValues=startValues, φ = λ[3])
            return mse(y[start:end], f[start:end])
        end
        lower_bound = zeros(3)
        upper_bound = ones(3)
        res = optimize(optPhi, lower_bound, upper_bound, [.8, .8, .8])
        λ = Optim.minimizer(res)
        f, D = holtLinear(y, λ[1], λ[2], startValues=startValues, h=h, φ = λ[3])
        return (forecast=f, model=D, parameters=λ)
    else
        function opt(λ)
            f, _ = holtLinear(y, λ[1], λ[2], startValues=startValues, φ=φ)
            return mse(y[start:end], f[start:end])
        end
        lower_bound = zeros(2)
        upper_bound = ones(2)
        res = optimize(opt, lower_bound, upper_bound, [.8, .8])
        λ = Optim.minimizer(res)
        f, D = holtLinear(y, λ[1], λ[2], startValues=startValues, h=h, φ = φ)
        return (forecast=f, model=D, parameters=λ)
    end 
end


"""
    holtWinters(y :: Vector, λ₁ :: Real, λ₂ :: Real, λ₃ :: Real, s::Int;
                     startValues::Tuple=(), h::Int = 0, model=:add, φ::Real = 1.0)

Holt-Winters method with smoothing factors λ₁ (level), λ₂ (slope) and λ₃ (season)
and damping factor φ (\\varphi) given vector y (N x 1) and h forecast periods.

The number of seasons s has to be specified (s>0), else use TrendDecomposition.holtLinear. 

Start values of level and slope can be given as tuple (l₀, b₀),
otherwise output values start at period 3.

Returns a tuple with (f, D), where f is a ((N+h) x 1) vector with the
forecast values and D is a (N x 3) matrix with the level, slope, season values as columns.

"""
function holtWinters(y :: Vector, λ₁ :: Real, λ₂ :: Real, λ₃ :: Real, s::Int;
                     startValues::Tuple=(), h::Int = 0, model=:add, φ::Real = 1.0)
    T = length(y)
    F = zeros(T + h)
    lvl = zeros(T)
    b = zeros(T)

    S = zeros(T)
    
    p = (mod(s, 2) == 0 ? s + 1 : s)
    S[1:s] = maSeason(y[1:(2*s)] .- rollingAverage(y[1:(2*s)], p), s)
    
    start :: Int = 3
    
    if startValues == ()
        b[1] = lvl[1] = NaN
        b[2] = y[2] - y[1]
        lvl[2] = y[2]
        F[1] = NaN
        F[2] = NaN
    else
        lvl[1] = λ₁*(y[1] - S[1]) + (1. - λ₁)*(startValues[1] + startValues[2] * φ)
        b[1] = λ₂ * (lvl[1] - startValues[1]) + (1. - λ₂) * startValues[2]  * φ
        S[1] = λ₃ * (y[1] - lvl[1]) + (1. - λ₃) * S[1]
        F[1] = startValues[1] + startValues[2]  * φ + S[1]
        start = 2
    end

    if model == :add
         S[1:s] = maSeason(y[1:(2*s)] .- rollingAverage(y[1:(2*s)], p), s)
        for t in start:T
            F[t] = lvl[t-1] + b[t-1] * φ + S[(t-s < 1 ? t : t-s)]
            lvl[t] = λ₁*(y[t] - S[(t-s < 1 ? t : t-s)]) + (1. - λ₁)*(lvl[t-1] + b[t-1] * φ)
            b[t] = λ₂ * (lvl[t] - lvl[t-1]) + (1. - λ₂) * b[t-1] * φ
            S[t] = λ₃ * (y[t] - lvl[t]) + (1. - λ₃) * S[(t-s < 1 ? t : t-s)]
        end

        if φ == 1.0
            for j in 1:h
                idx = (mod(j, s) == 0 ? s : mod(j, s))
                F[T + j] = lvl[T] + b[T] * j + S[(T-s+1):T][idx]
            end
        else
            for j in 1:h
                idx = (mod(j, s) == 0 ? s : mod(j, s))
                F[T + j] = lvl[T] + b[T] * sum([φ^i for i in 1:j]) + S[(T-s+1):T][idx]
            end
        end

        
    elseif model == :mul
        S[1:s] = maSeason(y[1:(2*s)] ./ rollingAverage(y[1:(2*s)], p), s)
        for t in start:T
            F[t] = (lvl[t-1] + b[t-1] * φ) * S[(t-s < 1 ? t : t-s)]
            lvl[t] = λ₁*(y[t] / S[(t-s < 1 ? t : t-s)]) + (1. - λ₁)*(lvl[t-1] + b[t-1] * φ)
            b[t] = λ₂ * (lvl[t] - lvl[t-1]) + (1. - λ₂) * b[t-1] * φ
            S[t] = λ₃ * (y[t] / lvl[t]) + (1. - λ₃) * S[(t-s < 1 ? t : t-s)]
        end

        if φ == 1.0
            for j in 1:h
                idx = (mod(j, s) == 0 ? s : mod(j, s))
                F[T + j] = (lvl[T] + b[T] * j) * S[(T-s+1):T][idx]
            end
        else
            for j in 1:h
                idx = (mod(j, s) == 0 ? s : mod(j, s))
                F[T + j] = (lvl[T] + b[T] * sum([φ^i for i in 1:j])) * S[(T-s+1):T][idx]
            end
        end

    end

    return (F, hcat(lvl, b, S))
end



mse(y, f) = sum((y .- f).^2) / length(y)


"""
    holtWinters(y :: Vector, s::Int;
                     startValues::Tuple=(), h::Int = 0, model=:add, φ::Union{Real, Nothing} = nothing)

Holt-Winters method with damping factor φ (\\varphi) given vector y (N x 1) and h forecast periods.

The number of seasons s has to be specified (s>0), else use TrendDecomposition.holtLinear. 

The smoothing parameters are estimated by minimizing the loss
function using Optim.jl. For φ equal nothing, its value is also determined by the optimization algorithm.

Start values of level and slope can be given as tuple (l₀, b₀),
otherwise output values start at period 3.

Returns a tuple with (f, D, p), where f is a ((N+h) x 1) vector with the
forecast values and D is a (N x 3) matrix with the level, slope, season values as columns,
and p is a vector containing the estimated parameter(s).
"""
function holtWinters(y :: Vector, s::Int;
                     startValues::Tuple=(), h::Int = 0, model=:add, φ::Union{Real, Nothing} = nothing)
    
    start :: Int = (startValues == () ? 3 : 1)
    if φ == nothing
        function optPhi(λ)
            f, _ = holtWinters(y, λ[1], λ[2], λ[3], s, startValues=startValues, model=model, φ = λ[4])
            return mse(y[start:end], f[start:end])
        end
        lower_bound = zeros(4)
        upper_bound = ones(4)
        res = optimize(optPhi, lower_bound, upper_bound, [.8, .8, .8, .8])
        λ = Optim.minimizer(res)
        f, D = holtWinters(y, λ[1], λ[2], λ[3], s, startValues=startValues, h=h, model=model, φ = λ[4])
        return (forecast=f, model=D, parameters=λ)
    else
        function opt(λ)
            f, _ = holtWinters(y, λ[1], λ[2], λ[3], s, startValues=startValues, model=model, φ=φ)
            return mse(y[start:end], f[start:end])
        end
        lower_bound = zeros(3)
        upper_bound = ones(3)
        res = optimize(opt, lower_bound, upper_bound, [.8, .8, .8])
        λ = Optim.minimizer(res)
        f, D = holtWinters(y, λ[1], λ[2], λ[3], s, startValues=startValues, h=h, model=model, φ = φ)
        return (forecast=f, model=D, parameters=λ)
    end 
end
