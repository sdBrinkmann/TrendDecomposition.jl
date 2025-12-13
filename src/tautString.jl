

"""
    tautStringFit(y :: Vector, C :: Real)

Computes the taut string for time series y where C determines the diameter of
the tube.

Returns (string, (x, y)) with string being the fitted taut string and
the x and y coordinates of the knots. 

"""
function tautStringFit(y :: Vector, C :: Real)

    n :: Int = length(y)
    y_integrated = vcat(0, cumsum(y) ./ n)
    
    lower = y_integrated .- (C / sqrt(n))
    upper = y_integrated .+ (C / sqrt(n))


    p1, k1 = tautString(lower, upper, C, lower[1], lower[n])

    #= Optimization not yet implemented
    p2, k2 = tautString(lower, upper, C, lower[1], upper[n])
    p3, k3 = tautString(lower, upper, C, upper[1], lower[n])
    p4, k4 = tautString(lower, upper, C, upper[1], upper[n])
    =#

    string = Vector{Float64}(undef, n)

    no_knots = length(p1)

    x_range = range(0, 1, (n+1))

 
    factor = 1 / n
    
      for i in 1:(no_knots-1)
        slope = (k1[i+1] - k1[i]) / ((p1[i+1] - p1[i]) * factor)
        string[p1[i]:(p1[i+1]-1)] .= slope

    end

    
    return (string, (p1, k1))
end


function tautString(lower :: Vector, upper :: Vector, C :: Real, startValue :: Real, endValue :: Real)     
    n :: Int = length(lower)
    
    #lower[1] = upper[1] = startValue
    #lower[n] = upper[n] = endValue

    u_knots = zeros(n)
    l_knots = zeros(n)
    s_knots = zeros(n)

    u_points = ones(Int, n)
    l_points = ones(Int, n)
    s_points = ones(Int, n)

    u_knots[1] = l_knots[1] = s_knots[1] = startValue   #lower[1]
    
    u_last :: Int = l_last :: Int = 1
    u_index :: Int = l_index :: Int = 1

    s_index :: Int = 1

    u_base :: Int = l_base :: Int = 0
    
    i :: Int = 2

    while i < n
        
        if u_index < i
            u_index +=1
            #u_last = pool!(>, u_knots, u_points,
            #               u_last, u_index, upper[u_index])
            u_last = pool!(>, view(u_knots,1+u_base:n), view(u_points, 1+u_base:n),
                           u_last, u_index, upper[u_index])
        end

        if l_index < i
            l_index += 1
            #l_last = pool!(<, l_knots, l_points,
            #               l_last, l_index, lower[l_index])
            l_last = pool!(<, view(l_knots, 1+l_base:n), view(l_points, 1+l_base:n),
                           l_last, l_index, lower[l_index])
        end

        # sxᵢ(0+) > svᵢ(0+)

        if ((u_knots[2+u_base] - u_knots[1+u_base]) / (u_points[2+u_base] - u_points[1+u_base])) >
            ((l_knots[2+l_base] - l_knots[1+l_base]) / (l_points[2+l_base] - l_points[1+l_base]))
            i += 1
        else
            s_index += 1
            ux :: Int = 2 + u_base
            lx :: Int = 2 + l_base
            #println("$i: false, $(u_points[ux]) $(l_points[lx])")
            if u_points[ux] < l_points[lx]
                s_knots[s_index] = u_knots[ux]
                s_points[s_index] = u_points[ux]
                i = u_points[ux] + 1
                l_last = 1
                l_index = u_points[ux]
                l_knots[1] = u_knots[ux]
                l_points[1] = u_points[ux]

                u_last -= 1
                
                u_base += 1
                l_base = 0
            else
                s_knots[s_index] = l_knots[lx]
                s_points[s_index] = l_points[lx]
                i = l_points[lx] + 1

                u_last = 1
                u_index = l_points[lx]
                u_knots[1] = l_knots[lx]
                u_points[1] = l_points[lx]

                l_last -= 1
                
                l_base += 1
                u_base = 0
                
            end
            
        end
    end

    s_index += 1
    s_points[s_index] = n
    s_knots[s_index] = lower[n]
    return (s_points[1:s_index], s_knots[1:s_index])
end



## Equidistant pool

function pool!(compare :: Function, knots :: AbstractVector, points :: AbstractVector, last :: Int, x :: Real, y :: Real)
    if last == 1
        knots[2] = y
        points[2] = x
        return 2
    else
        slope1 = (y - knots[last]) / (x - points[last])
        slope2 = (knots[last] - knots[last - 1]) / (points[last] - points[last-1])

        if compare(slope1, slope2)
            last += 1
            knots[last] = y
            points[last] = x
            return last
        else
            pool!(compare, knots, points, last - 1, x, y)
        end
    end
end


"""
    leastConcaveMajorant(y :: Vector)

Computes least concave majorant (lcm) of series y using the
pool-adjecent-violators algorithm

Returns the coordinates of the knots as a tuple (x, y) 
"""
function leastConcaveMajorant(y :: Vector)
    n = length(y)

    knots = zeros(n)
    points = ones(Int, n)
    last = 1

    knots[1] = y[1]
    
    for i in 2:n
        last = pool!(<, knots, points, last, i, y[i])
    end

    return (points[1:last], knots[1:last])
end


"""
    greatestConvexMinorant(y :: Vector)

Computes greatest convex minorant of series y using the
pool-adjecent-violators algorithm

Returns the coordinates of the knots as a tuple (x, y) 
"""
function greatestConvexMinorant(y :: Vector)
    n = length(y)

    knots = zeros(n)
    points = ones(Int, n)
    last = 1

    knots[1] = y[1]
    
    for i in 2:n
        last = pool!(>, knots, points, last, i, y[i])
    end

    return (points[1:last], knots[1:last])
end



