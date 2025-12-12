




function softThreshold(a, κ)
    if a > κ
        return (a - κ)
    elseif a < -κ
        return (a + κ)
    else
        return 0
    end
end



"""
    trendADMM(y :: Vector, λ :: Real; m = 2, max_iter = 100, criterion = true)

Trend filtering of time series data y by using the L1 penalty with regularization parameter λ
and using the m'th difference.
The estimation procedure uses the Alternating Direction Method of Multipliers (ADMM) to
reach a numerical solution.

If criterion = true the function stop when reaching a criterion's stopping condition,
otherwise algorithm iterates max_iter times.

The function returns the estimated trend component.

"""
function trendADMM(y :: Vector, λ :: Real; m = 2, max_iter = 100, criterion = true)
    n = length(y) 
    z = zeros(n)
    u = zeros(n)
    τ = zeros(n)

    if λ < 0
        throw(DomainError(λ, "Only positive values can be used for λ"))
    elseif max_iter < 1
        throw(DomainError(max_iter, "max_iter must be positive"))        
    end

    
    
    stopped = false
    
    ϵ_abs = 1.e-4
    ϵ_rel = 1.e-2
    
    D = difference_matrix(n, m, full=true)
    ρ = λ

    DD = D * D'

    for i in 1:max_iter
        τ = (I + ρ * DD) \ (y + ρ * D * (z + u))

        z_prv = z
        z = softThreshold.(D'*τ - u, λ / ρ)
        #println(D'*τ - u)
        
        u = u + z - D' * τ

        # Test Stop criterion
        r_norm = sqrt(sum((z - D'* τ).^2))
        #println("z: $z,    z_old: $z_prv")
        s_norm = sqrt(sum((ρ * D' * (z - z_prv)).^2))

        ϵ_primal = sqrt(n) * ϵ_abs + ϵ_rel * max(norm(D'*τ), norm(z))
        ϵ_dual = sqrt(n) * ϵ_abs + ϵ_rel * norm(ρ * D * u)

        #println("$i    $r_norm    $ϵ_primal    $s_norm    $ϵ_dual")
        
        if criterion && (r_norm < ϵ_primal  && s_norm < ϵ_dual)
            println("Stopping criterion met at iteration $i")
            stopped = true
            break
        end       
    end

    if stopped == false
        println("Stop criterion was not reached, stopped after max iterations $(max_iter)")
    end
    
    return τ
end


"""
    trendL1Filter()

Placeholder for trendL1Filter extension, when using TrendDecomposition together with other Julia packages
like Convex.jl and SCS.jl.

trendL1Filter serves as a wrapper function to enable the usage of serveral different methods in estimating
L1 trend filter. 
"""
function trendL1Filter() end

    

