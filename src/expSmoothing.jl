



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




function brownLinear(y :: Vector, λ :: Real; h::Int = 0)
    T = length(y)
    F = zeros(T)
    lvl = zeros(T)
    b = zeros(T)

    #if startValue == false
        b[2] = y[2] - y[1]
        lvl[2] = y[2]
    #end

    for t in 3:T
        F[t] = lvl[t-1] + b[t-1]
        lvl[t] = lvl[t-1] + b[t-1] + (1. - λ^2) * (y[t] - F[t])
        b[t] = b[t-1] + (1. - λ)^2 * (y[t] - F[t])
    end

    for j in 1:h
        F[T + j] = lvl[T] + b[T] * j
    end

    return (F, hcat(lvl, b))
end





function holtLinear(y :: Vector, λ₁ :: Real, λ₂ :: Real;
                    startValue::Bool=false, h::Int = 0)
    T = length(y)
    F = zeros(T + h)
    lvl = zeros(T)
    b = zeros(T)

    if startValue == false
        b[2] = y[2] - y[1]
        lvl[2] = y[2]
    end

    for t in 3:T
        lvl[t] = λ₁*y[t] + (1- λ₁)*(lvl[t-1] + b[t-1])
        b[t] = λ₂ * (lvl[t] - lvl[t-1]) + (1 - λ₂) * b[t-1]
        F[t] = lvl[t-1] + b[t-1]
    end

    for j in 1:h
        F[T + j] = lvl[T] + b[T] * j
    end

    return (F, hcat(lvl, b))
end





function holtWinters(y :: Vector, λ₁ :: Real, λ₂ :: Real, λ₃ :: Real, s::Int;
                     startValue::Bool=false, h::Int = 0, model=:add)
    T = length(y)
    F = zeros(T + h)
    lvl = zeros(T)
    b = zeros(T)

    S = zeros(T)
    
    p = (mod(s, 2) == 0 ? s + 1 : s)
    S[1:s] = maSeason(y[1:(2*s)] .- rollingAverage(y[1:(2*s)], p), s)
    
    if startValue == false
        b[2] = y[2] - y[1]
        lvl[2] = y[2]
    end

    if model == :add
         S[1:s] = maSeason(y[1:(2*s)] .- rollingAverage(y[1:(2*s)], p), s)
        for t in 3:T
            lvl[t] = λ₁*(y[t] - S[(t-s < 1 ? t : t-s)]) + (1. - λ₁)*(lvl[t-1] + b[t-1])
            b[t] = λ₂ * (lvl[t] - lvl[t-1]) + (1. - λ₂) * b[t-1]
            S[t] = λ₃ * (y[t] - lvl[t]) + (1. - λ₃) * S[(t-s < 1 ? t : t-s)]
            F[t] = lvl[t-1] + b[t-1] + S[(t-s-1 < 1 ? t-1 : t-s-1)]
        end

        for j in 1:h
            idx = (mod(j, s) == 0 ? s : mod(j, s))
            F[T + j] = lvl[T] + b[T] * j + S[(T-s+1):T][idx]

        end
    elseif model == :mul
        S[1:s] = maSeason(y[1:(2*s)] ./ rollingAverage(y[1:(2*s)], p), s)
        for t in 3:T
            lvl[t] = λ₁*(y[t] / S[(t-s < 1 ? t : t-s)]) + (1. - λ₁)*(lvl[t-1] + b[t-1])
            b[t] = λ₂ * (lvl[t] - lvl[t-1]) + (1. - λ₂) * b[t-1]
            S[t] = λ₃ * (y[t] / lvl[t]) + (1. - λ₃) * S[(t-s < 1 ? t : t-s)]
            F[t] = (lvl[t-1] + b[t-1]) * S[(t-s-1 < 1 ? t-1 : t-s-1)]
        end

        for j in 1:h
            idx = (mod(j, s) == 0 ? s : mod(j, s))
            F[T + j] = (lvl[T] + b[T] * j) * S[(T-s+1):T][idx]
        end

    end

    return (F, hcat(lvl, b, S))
end

