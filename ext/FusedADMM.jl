
module FusedADMM

using TrendDecomposition, Lasso, LinearAlgebra



"""
    fusedADMM(y :: Vector, λ :: Real; m::Int = 2, max_iter::Int = 1000, ρ::Real=λ,
                ϵ_abs::Float64=1.e-4, ϵ_rel::Float64=1.e-2)

Trend filtering of time series data y by using the L1 penalty with regularization parameter λ
and using the m'th difference.
The estimation procedure uses the Alternating Direction Method of Multipliers (ADMM) in combination
with the fused lasso algorithm from the Lasso.jl package to reach a numerical solution. For m = 1 only the taut string algorithm
is used as an edge case solution. 


The function returns the estimated trend component.
"""
function TrendDecomposition.fusedADMM(y :: Vector, λ :: Real;
                                      m::Int = 2, max_iter::Int = 1000, ρ::Real=λ,
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
        string = coef(fit(FusedLasso, y, λ))
        return string
    else
        m -= 1
    end
    
    #ϵ_abs = 1.e-4
    #ϵ_rel = 1.e-2
    
    D = TrendDecomposition.difference_matrix(n, m, full=true)
    #ρ = λ

    DD = D * D'

    for i in 1:max_iter
        τ = (I + ρ * DD) \ (y + ρ * D * (z + u))

        z_prv = z
        z = coef(fit(FusedLasso, D'*τ - u, λ / ρ))
        
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
            println("stop criterion met at $i")
            break
        end
      
    end
    return τ
end




end #end module 
