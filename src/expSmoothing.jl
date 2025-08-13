



"""
    expSmoothing(y :: Vector, λ :: Real; startValue::Bool = true, v_0::Real = 0.0)

    Simple exponential smoothing with smoothing factor λ. If startValue equals true,
    v_0 will be used as initial value, otherwise the first element of y (y[1]) will be selected as
    starting value. Using y[1] as starting value will result in a different mathematical function.   
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
