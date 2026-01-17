




function softThreshold(a :: Float64, κ :: Float64)
    if a > κ
        return (a - κ)
    elseif a < -κ
        return (a + κ)
    else
        return 0
    end
end



"""
    trendADMM(y :: Vector, λ :: Real; m::Int = 2, max_iter::Int = 2000, ρ::Real=λ,
                ϵ_abs::Float64=1.e-4, ϵ_rel::Float64=1.e-3)

Trend filtering of time series data y by using the L1 penalty with regularization parameter λ
and using the m'th difference.
The estimation procedure uses the Alternating Direction Method of Multipliers (ADMM) to
reach a numerical solution.

The function returns the estimated trend component.

"""
function trendADMM(y :: Vector, λ :: Real;
                   m::Int = 2, max_iter::Int = 2000, ρ::Real=λ,
                   ϵ_abs::Float64=1.e-4, ϵ_rel::Float64=1.e-3)
    n = length(y) 
    z = zeros(Float64, n)
    u = zeros(Float64, n)
    τ = zeros(Float64, n)

    if λ < 0
        throw(DomainError(λ, "Only positive values can be used for λ"))
    elseif max_iter < 1
        throw(DomainError(max_iter, "max_iter must be positive"))
    elseif m < 1
        throw(DomainError(m, "Order of differentation must be greater than 0"))
    elseif ρ < 1.
        throw(DomainError(ρ, "Only values greater equal 1 can be used for ρ"))
    end

    ρ = Float64(ρ)
    
    
    #ϵ_abs = 1.e-4
    #ϵ_rel = 1.e-3
    
    D = difference_matrix(n, m, full=true)
    #ρ = λ

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
     
        if r_norm < ϵ_primal  && s_norm < ϵ_dual
            println("Stop criterion met at iteration $i")
            return τ
        end
        
    end

    println("Stop criterion was not reached, stopped after max iterations $(max_iter)")

    
    return τ
end


"""
    tautADMM(y :: Vector, λ :: Real;
                m::Int = 2, max_iter::Int = 2000, ρ::Real=λ, opt::Bool = false,
                ϵ_abs::Float64=1.e-4, ϵ_rel::Float64=1.e-2)

Trend filtering of time series data y by using the L1 penalty with regularization parameter λ
and using the m'th difference.
The estimation procedure uses the Alternating Direction Method of Multipliers (ADMM) in combination
with the taut string algorithm to reach a numerical solution. For m = 1 only the taut string algorithm
is used as an edge case solution. 


The function returns the estimated trend component.
"""
function tautADMM(y :: Vector, λ :: Real;
                  m::Int = 2, max_iter::Int = 2000, ρ::Real=λ, opt::Bool = false,
                  ϵ_abs::Float64=1.e-4, ϵ_rel::Float64=1.e-2)
    n = length(y) 
    z = zeros(Float64, n)
    u = zeros(Float64, n)
    τ = zeros(Float64, n)

    if λ < 0
        throw(DomainError(λ, "Only positive values can be used for λ"))
    elseif max_iter < 1
        throw(DomainError(max_iter, "max_iter must be positive"))
    elseif m < 1
        throw(DomainError(m, "Order of differentation must be greater than 0"))
    elseif ρ < 1.
        throw(DomainError(ρ, "Only values greater equal 1 can be used for ρ"))
    end
    
    if m == 1
        string, _ = tautStringFit(y, λ, optimize=opt)
        return string
    else
        m -= 1
    end

    ρ = Float64(ρ)
    
    #ϵ_abs = 1.e-4
    #ϵ_rel = 1.e-2
    
    D = difference_matrix(n, m, full=true)

    DD = D * D'

    for i in 1:max_iter
        τ = (I + ρ * DD) \ (y + ρ * D * (z + u))

        z_prv = z
        z, _ = tautStringFit(D'*τ - u, λ / ρ, optimize=opt)
        #println(D'*τ - u)
        
        u = u + z - D' * τ

        # Test Stop criterion
        r_norm = sqrt(sum((z - D'* τ).^2))
        #println("z: $z,    z_old: $z_prv")
        s_norm = sqrt(sum((ρ * D' * (z - z_prv)).^2))

        ϵ_primal = sqrt(n) * ϵ_abs + ϵ_rel * max(norm(D'*τ), norm(z))
        ϵ_dual = sqrt(n) * ϵ_abs + ϵ_rel * norm(ρ * D * u)

        #println("$i    $r_norm    $ϵ_primal    $s_norm    $ϵ_dual")
        
        if (r_norm < ϵ_primal  && s_norm < ϵ_dual)
            println("Stop criterion met at $i")
            return τ
        end
        
    end

    println("Stop criterion was not reached, stopped after max iterations $(max_iter)")
    
    return τ
end




"""
    trendL1Filter(y :: Vector, λ :: Real; m = 2, max_iter=20, method = :ADMM)

Placeholder for trendL1Filter extension, when using TrendDecomposition together with other Julia packages
like Convex.jl and SCS.jl.

This function provides the generic use of serveral optimization methods to compute a numerical solution.
Following methods are implmented:
:ADMM -> alternating direction method of multipliers
:ConvexSCS -> SCS solver with Convex.jl. Prerequisite! Import necessary modules with: using Convex, SCS

The function returns the estimated trend component
"""
function trendL1Filter() end



    

"""
    fusedADMM(y :: Vector, λ :: Real; m::Int = 2, max_iter::Int = 1000, ρ::Real=λ,
                ϵ_abs::Float64=1.e-4, ϵ_rel::Float64=1.e-2)

Placeholder for FusedADMM extension, when using TrendDecomposition together with Lasso.jl package.

Trend filtering of time series data y by using the L1 penalty with regularization parameter λ
and using the m'th difference.
The estimation procedure uses the Alternating Direction Method of Multipliers (ADMM) in combination
with the fused lasso algorithm from the Lasso.jl package to reach a numerical solution. For m = 1 only the taut string algorithm
is used as an edge case solution. 


The function returns the estimated trend component.
"""
function fusedADMM() end
