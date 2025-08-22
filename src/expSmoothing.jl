



"""
    expSmoothing(y :: Vector, λ :: Real; startValue::Bool = true, v_0::Real = 0.0)

Simple exponential smoothing with smoothing factor λ. If startValue equals true,
v_0 will be used as initial value, otherwise the first element of y (y[1])
will be selected as starting value.
Using y[1] as starting value will result in a different mathematical function.   
"""
function expSmoothing(y :: Vector, λ :: Real; start_value::Bool = true, v_0::Real = 0.0)
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
φ (\varphi) given vector y (N x 1) and h forecast periods.

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
    holtLinear(y :: Vector, λ₁ :: Real, λ₂ :: Real;
                    startValues::Tuple=(), h::Int = 0, φ::Real = 1.0)

Holt's linear trend method with smoothing factors λ₁ (level) and λ₂ (slope)
and damping factor φ (\varphi) given vector y (N x 1) and h forecast periods.

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
    holtWinters(y :: Vector, λ₁ :: Real, λ₂ :: Real, λ₃ :: Real, s::Int;
                     startValues::Tuple=(), h::Int = 0, model=:add, φ::Real = 1.0)

Holt-Winters method with smoothing factors λ₁ (level), λ₂ (slope) and λ₃ (season)
and damping factor φ (\varphi) given vector y (N x 1) and h forecast periods.

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
        F[1:2] .= NaN
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
            lvl[t] = λ₁*(y[t] - S[(t-s < 1 ? t : t-s)]) + (1. - λ₁)*(lvl[t-1] + b[t-1] * φ)
            b[t] = λ₂ * (lvl[t] - lvl[t-1]) + (1. - λ₂) * b[t-1] * φ
            S[t] = λ₃ * (y[t] - lvl[t]) + (1. - λ₃) * S[(t-s < 1 ? t : t-s)]
            F[t] = lvl[t-1] + b[t-1] * φ + S[(t-s-1 < 1 ? t-1 : t-s-1)]
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
            lvl[t] = λ₁*(y[t] / S[(t-s < 1 ? t : t-s)]) + (1. - λ₁)*(lvl[t-1] + b[t-1] * φ)
            b[t] = λ₂ * (lvl[t] - lvl[t-1]) + (1. - λ₂) * b[t-1] * φ
            S[t] = λ₃ * (y[t] / lvl[t]) + (1. - λ₃) * S[(t-s < 1 ? t : t-s)]
            F[t] = (lvl[t-1] + b[t-1] * φ) * S[(t-s-1 < 1 ? t-1 : t-s-1)]
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


