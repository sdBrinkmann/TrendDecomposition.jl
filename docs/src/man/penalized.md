# Penalized Smoothing

The term penalized smoothing here denotes a class of smoothers  that minimize a criterion 
subject to a penalizing term. The methods
introduced here are in fact one of the earliest uses in statistics of a 
penalizing term first introduced in 1898 by George Bohlmann[^Bohl99] and it is
the first smoothing procedure posed as an optimization problem. 

[^Bohl99]:
	> Bohlmann, G. (1899). Ein Ausgleichungsproblem. Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse, 1899, 260-271.

Given our time series ``y``, it is assumed that it can be decomposed into
```math
	y_t = \tau_t + c_t,
```
where ``\tau_t`` is the trend component and ``c_t`` is the cycle component at time t. 
The criterion we want to minimize is the squared difference between the 
trend component and the corresponding data point ``(y_t - \tau_t)^2``.

For the penalizing term we differentiate between the traditionally used ``l_2``
norm and more novel approaches using the ``l_1`` norm, similar to L1 and L2 regularization in regression analysis. 


## Whittaker-Henderson Smoothing

Using the squared m'th difference of ``\tau_t``, written here as ``(\Delta^{m} \tau_t)^2``,
as penalizing term and ``\lambda > 0`` as regularization parameter, the trend
component can be estimated as solution to the optimization problem:

```math
\^\tau_t = 	\operatorname*{arg\,min}_\tau \Bigg\{ \sum_{t=1}^{n} (y_t - \tau_t)^2 + \lambda \sum_{t=m+1}^{n} 
(\Delta^m \tau_t)^2\Bigg\}.
```

All functions below, `bohlmannFilter`, `hpFilter`, `bhpFilter`, use the algebraic solution to the minimization problem. 
The ``\lambda`` value has to be selected beforehand by the user.
For ``\lambda \to 0`` the solution converges to the original data and for ``\lambda \to \infty`` the trend becomes a straight line. 


```@docs
bohlmannFilter
```

## Hodrick-Prescott (HP) Filter
The case using second difference with m = 2 was proposed by Leser (1961)[^Leser61] for trend estimation and the above equation becomes 

[^Leser61]:
	> Leser, C. E. V. (1961). A Simple Method of Trend Construction. Journal of the Royal Statistical Society. Series B (Methodological), 23(1), 91–107.

```math
\^\tau_t = 	\operatorname*{arg\,min}_\tau \Bigg\{ \sum_{t=1}^{n} (y_t - \tau_t)^2 + \lambda \sum_{t=m+1}^{n} 
(\tau_{t+1} - 2 \tau_t + \tau_{t-1})^2\Bigg\}.
```

In the literature this procedure is also known as the Hodrick-Prescott filter.  In general, the following parameters are recommended for ``\lambda``

| Basis     | Period length     | λ           |
| ------    | -----------       | ----------  |
| Annual    | 1                 | 100         |
| Quarterly | 4                 | 1600        |
| Monthly   | 12                | 14,600      |





```@docs
hpFilter(x::Vector, λ::Real)
```

## Boosted HP Filter

In order to provide a better trend estimation, according to Philips and Shi (2021)[^PhilShi21]
the HP Filter can be applied repeatedly over the cycle estimate(s). This idea
is akin to ``L_2``-boosting in machine learning, hence the method is called boosted HP
(bHP) filter. 

[^PhilShi21]:
	> Phillips, P.C.B. and Shi, Z. (2021), BOOSTING: WHY YOU CAN USE THE HP FILTER. International Economic Review, 62: 521-570.

`hpFilter(x::Vector, λ::Real, iter::Int)` allows the repeated application of the HP filter by a fixed amount of times determind by the parameter `iter`.

```@docs
hpFilter(x::Vector, λ::Real, iter::Int)
```

`bhpFilter` implements an iterative approach, where a stopping criterion determines
the amout of iterations.

```@docs
bhpFilter
```

# ``l_1`` Trend Filtering

## Taut String - Piecewise constant 

The taut string algorithm by Davies and Kovac (2001)[^DaKo01] can be seen as an efficient O(n) solution to the following minimization problem:

[^DaKo01]:
	> Davies, P. L., & Kovac, A. (2001). Local extremes, runs, strings and multiresolution. The Annals of Statistics, 29(1), 1-65

```math
\^\tau_t = 	\operatorname*{arg\,min}_\tau \Bigg\{ \sum_{t=1}^{n} (y_t - \tau_t)^2 + \lambda \sum_{t=m+1}^{n} |\tau_t - \tau_{t-1}|\Bigg\},
```

where ``\tau`` is again the trend component. The result is a piecewise constant function and it
can be shown that is uses the smallest number of knots.

Alternatively, using the term total variation (TV) of a function f the taut string approach minimizes
following function:

```math
 \sum_{t=1}^{n} (y_t - f(x_t))^2 + \lambda TV(f),
```
also known as total variation regularization or denoising.

Instead of the parameter ``\lambda``, the parameter ``C`` is used to determine the radius of the tube
in the taut string method, but in practise ``C`` is equivalent to the parameter ``\lambda`` in the
alternative estimation methods below (ADMM and fused Lasso).

```@docs
tautStringFit
```

This function implements a suggestion by the authors of the algorithm to vary the starting and end points in order to find the best fit. With `optimize = true` various end point are tested and the best outcome is choosen. This is a pure experimental feature. Run julia with `julia --threads 5` to take advantage of multi-threading. 

## Alternative Direction Method of Multipliers (ADMM)

To use ADMM as a solution method to total variation minimization was proposed by Boyd et al. (2011)[^Boyd11]

[^Boyd11]:
	> Boyd, S., Neal, P., Eric, C., Borja, P., & Jonathan, E. (2011). Distributed optimization and statistical learning via the alternating direction method of multipliers. Foundations and Trends® in Machine learning, 3(1), 1-122.

and is here implemented as a general solver for any order of l1 trend filtering. 

```@docs 
trendADMM
```

### Alternatives to Soft-Thresholding

Ramas and Tibshirani (2016)[^RaTi12] proposed using the fused Lasso instead of soft-thresholding to archieve faster convergence and a more numerical robust solution. A second alternative is offered here using the taut string algorithm instead of the fused Lasso.

[^RaTi12]:
	> Ramdas, A., & Tibshirani, R. J. (2016). Fast and flexible ADMM algorithms for trend filtering. Journal of Computational and Graphical Statistics, 25(3), 839-858


#### Taut String

```@docs
tautADMM
```

In the same vein as with `tautStringFit` above, here `opt=true` can be used to experiment with various 
end points in the taut string algorithm.

#### Fused Lasso
This function is offered as an extension for the package Lasso.jl. Please use with `using Lasso`.
```@docs
fusedADMM
```


### Note on augmented Lagrangian parameter ρ

Using ADMM with soft-thresholding or taut string the augmented Lagrangian parameter ``\rho`` (written
as `\rho<TAB>`) can be tuned to archieve a faster convergence in practise. As a parameter for the subroutine using either soft-thresholding or taut string ``\lambda / \rho`` is used with ``\rho = \lambda`` as default (same default for fused Lasso). Using a ``\rho`` greater than or rather a multiple of ``\lambda`` can increase the convergence rate significantly. 
