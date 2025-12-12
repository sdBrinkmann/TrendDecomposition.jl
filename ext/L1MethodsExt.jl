
module L1MethodsExt


using TrendDecomposition, Convex, SCS



function trendConvexJL(y :: Vector, λ :: Real; m=2)
    n = length(y)
    D = TrendDecomposition.difference_matrix(n, m)
    τ = Variable(n)
    ## problem = minimize(0.5 * v'*D'*D*v-y'*D*v, [v >= -λ, v <= λ])
    problem = minimize(0.5*sumsquares(y-τ)+ λ*norm(D'*τ,1))
    solve!(problem, SCS.Optimizer; silent = true)
    return τ.value
end

"""
    trendL1Filter(y :: Vector, λ :: Real; m = 2, max_iter=20, method = :ADMM)

General trend filtering of time series data y by using the L1 penalty with regularization
parameter λ in application of Bohlmann Filter or Whittaker-Henderson smoothing using the m'th difference.

This function provides the generic use of serveral optimization methods to compute a numerical solution.
Following methods are implmented:
:ADMM -> alternating direction method of multipliers
:ConvexJL -> SCS solver with Convex.jl. Prerequisite! Import necessary modules with: using Convex, SCS

The function returns the estimated trend component
"""
function TrendDecomposition.trendL1Filter(y :: Vector, λ :: Real; m = 2, max_iter=100, method = :ADMM)

    if method == :ADMM
        return trendADMM(y, λ, m=m, max_iter=max_iter)
    elseif method == :ConvexJL
        return trendConvexJL(y, λ, m=m)
    else
        println("Invalid method, valid choices as :ADMM or :ConvexJL")
        return []
    end

end



end #end module
