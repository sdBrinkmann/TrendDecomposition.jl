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
