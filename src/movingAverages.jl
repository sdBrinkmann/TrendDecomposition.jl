

import Base.sum
"""
    rollingAverage(y :: Vector, p :: Int; centered::Bool = true, discard::Bool = false, offset::Int = 0)

Computes the moving average of the vector y with p data points.

With centered equal true, p has to be an odd number so that an equal number of leads and lags can be used. By default
when centered = false, the computation includes the datum plus its most recent (p-1) lagged value.
This behavior can be changed by including an offset; a positive offset includes one additonal lead at
the expense of the oldest lag value

For discard=false the function rolls over all data, even if not all p data points can be used.
   
"""
function rollingAverage(y :: Vector, p :: Int; centered::Bool = true, discard::Bool = false, offset::Int = 0)
    len = length(y) 
    @assert len > p
    ave = zeros(len)
    if centered == true
        @assert mod(p, 2) == 1
        half = div(p, 2)
        @assert abs(offset) <= half
        for t in 1:len
            index = max(1, t-half+offset):min(t+half+offset, len)
            if discard == true && length(index) != p
                ave[t] = NaN
            else
                ave[t] = sum(y[index]) / length(index)
            end
        end
    else
        @assert offset >= 0
        @assert offset < p
        for t in 1:len
            index = max(1, t-p+1+offset):min(t+offset, len)
            if discard == true && length(index) != p
                ave[t] = NaN
            else
                ave[t] = sum(y[index]) / length(index)
            end
        end
    end
    #=
    if discard == true
        if centered == true
            return ave[(half+1-offset):(end-half-offset)]
        else
            return ave[p-offset:end-offset]
        end
    else
        return ave
    end
    =#
    return ave
end


function sum(A :: Vector, weights :: Vector) 
    acc = 0
    for (a, w) in zip(A, weights)
        acc += a * w
    end
    return acc
end

"""
    rollingAverage(y :: Vector, p :: Int, weights :: Vector;
        centered::Bool = true, discard::Bool = false, offset::Int = 0)

Computes the moving average of the vector y with p data points weighted by vector w.

With centered equal true, p has to be an odd number so that an equal number of leads and lags can be used. By default
when centered = false, the computation includes the datum plus its most recent (p-1) lagged value.
This behavior can be changed by including an offset; a positive offset includes one additonal lead at
the expense of the oldest lag value.

For discard=false the function rolls over all data, even if not all p data points can be used,
in any case all used weights will sum to 1.
   
"""
function rollingAverage(y :: Vector, p :: Int, weights :: Vector;
                        centered::Bool = true, discard::Bool = false, offset::Int = 0)
    @assert sum(weights) ≈ 1
    len_y = length(y)
    len_w = length(weights)
    @assert len_w == p   
    ave = zeros(length(y))
    if centered == true
        @assert mod(p, 2) == 1
        half = div(p, 2)
        @assert abs(offset) <= half
        center = half + 1 + (offset * (-1))
        for t in 1:len_y
            first = max(1, t-half+offset)
            last = min(len_y, t+half+offset)
            index = first:last
            if (length(index) == p)
                ave[t] = sum(y[index], weights)
            else
                ave[t] = sum(y[index],
                             weights[(center + first - t):(center + last - t)]
                             ./ sum(weights[(center + first - t):(center + last - t)]))
            end
        end
    else
        @assert offset >= 0
        @assert offset < p
        for t in 1:len_y
            first = max(1, t-p+1+offset)
            last = min(t+offset, len_y)
            index = first:last
            if length(index) == p
                ave[t] = sum(y[index], weights)
            else
                if discard == true
                    ave[t] = NaN
                else
                    shift = offset * (-1)
                    ave[t] = sum(y[index],
                                 weights[(first-t+p+shift):(last-t+p+shift)] ./
                                     sum(weights[(first-t+p+shift):(last-t+p+shift)]))
                end
                #=
                ave[t] =
                sum(y[index],
                weights[(len_w - length(index)+1):(end)] / sum(
                weights[(len_w - length(index)+1):(end)]))
                =#
            end
        end
    end
    #=
    if discard == true
        if centered == true
            return ave[(half+1-offset):(end-half-offset)]
        else
            return ave[p-offset:end-offset]
        end
    else
        return ave
    end
    =#
    return ave
end

"""
    maSeason(y :: Vector, seasons :: Int; repeating::Bool = false)

Given the number of seasons, the function computes the average value of each
season component. This method works better for detrended data.

With reapeating equal true the function will repeat the results until the output vector
has the same length as the input vector y.
"""
function maSeason(y :: Vector, seasons :: Int; repeating::Bool = false)
    len = length(y)
    S = zeros(seasons)
    j = 0
    for t in 1:len
        j += 1
        S[j] += y[t]

        if mod(j, seasons) == 0
            j = 0
        end
    end
              
    if mod(len, seasons) == 0
        S = S ./ (len / seasons)
    else
        for i in 1:mod(len, seasons)
            S[i] = S[i] / (div(len, seasons) + 1)
        end
        for i in (mod(len, seasons)+1):seasons
            S[i] = S[i] / div(len, seasons)
        end
    end

    if repeating == false
        return S
    else
        if mod(len, seasons) == 0
            return repeat(S, div(len, seasons))
        else
            return repeat(S, div(len, seasons) + 1)[1:len]
            
        end
    end
end

"""
    maDecompose(y :: Vector, seasons :: Int; combine::Bool = false)

Decomposes the time series y into trend, season and the remaining noise components.

Either assumes an additive model (:add) or a multiplicative model (:mul).

Replicates the R decompose{stats} function. For a more generic function implementation see
decompose of this package (TrendDecomposition.decompose).
"""
function maDecompose(y :: Vector, seasons :: Int; combine::Bool = false, model=:add)
    trend = rollingAverage(y, 5)
    if model == :add
        sean = maSeason(y - trend, seasons, repeating=true)
        idiot = y - trend - sean
    else
        sean = maSeason(y ./ trend, seasons, repeating=true)
        idiot = y ./ trend ./ sean
    end
    if combine == true
        return hcat(y, trend, sean, idiot)
    else
        return hcat(trend, sean, idiot)
    end
    #return idiot
end


