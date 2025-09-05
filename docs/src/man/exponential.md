# Exponential Smoothing

Exponential smoothing, also known as exponentially weighted moving average (EWMA),
allows us to model a time series as an additive or multiplicative composition,
in order to conduct forecasts. Here in its additive form it can be written as 

```math
	y_t = l_t + S_t + u_t,
```
where ``l`` denotes the level or trend component and ``S`` the seaonal component, the 
latter component is only included in the Holt-Winters method. 

The residual or noise term ``u_t`` is implicitly assumed and the parameters of the 
models down below are estimated by minimizing the squared residual terms, ``\^u_t = y_t - \^y_t``, 
summed over the entire time series.
The model makes otherwise no assumptions and is stated in its classic recursive fashion.

## Simple exponential smoothing
	
Simple exponential smoothing is a special case of weighted moving averages, where
weights decline exponentially. Implemented as a recursion it is know as 
exponentially weighted moving average (EWMA) defined as:


```math
    l_t = (1-\lambda) * l_{t-1} + \lambda * y_t
```
For the recursion to work a start value has to be defined. Given that the time series
``y_t`` is recorded for ``t = 1,...,T``, a start value ``y_0`` has to be selected.


```@docs
expSmoothing
```

## Holt's procedure

In order to account for a trend and to make forecasts a local slope ``b_t`` was introduced by Holt (1957)[^Holt57] , which is updated each period t. 
As for most time series it is unrealistic to assume that a linear line is to presist much long
into the future, a damping factor ``\varphi`` (\varphi) can be used and is set between ``0 < \varphi < 1``.

[^Holt57]:
	> Holt, Charles C. (1957). "Forecasting Trends and Seasonal by Exponentially Weighted Averages". Office of Naval Research Memorandum

```math
	\begin{aligned}
    l_t &= \lambda_1 * y_t + (1-\lambda_1) * (l_{t-1} + \varphi b_{t-1}) \\
	b_t &= \lambda_2 * (l_{t} + l_{t-1}) + (1-\lambda_2) * \varphi b_{t-1} \\
	\^y_t &= l_{t-1} + \varphi b_{t -1}
	\end{aligned}
```

```@docs
holtLinear(y :: Vector, λ₁ :: Real, λ₂ :: Real)
```


## Double exponential smoothing

By using the above introduced residual or error term ``\^u_t``, the level ``l_t`` and
slope ``b_t`` parameters can be estimated as:

```math
	\begin{aligned}
	l_t &= l_{t-1} + \varphi b_{t-1} + (1 - \lambda^2) \^u_t \\
	b_t &= \varphi b_{t-1} + (1 - \lambda)^{2} \^u_t  \\
	\^y_t &= l_{t-1} + \varphi b_{t -1}
	\end{aligned}
```
This recursion is known as double exponential smoothing and also named after R.G. Brown and
can been seen as a special case of Holt's procedure by setting ``\lambda_0 = 1 - \lambda^2``
and ``\lambda_1 = \frac{1 - \lambda}{1 + \lambda}``.


```@docs
brownLinear(y :: Vector, λ :: Real)
```


## Holt-Winters forecasting

In addition a seasonal component ``S_t`` can be modeled for number of seasons s.
This introduces an additional equation for ``S_t``: 

```math
	\begin{aligned}
    l_t &= \lambda_1 * (y_t - S_{t-s}) + (1-\lambda_1) * (l_{t-1} + \varphi b_{t-1}) \\
	b_t &= \lambda_2 * (l_{t} + l_{t-1}) + (1-\lambda_2) * b_{t-1} \varphi \\
    S_t &= \lambda_3 * (y_t + l_t) + (1-\lambda_3) * S_{t-s} \\
	\^y_{T+h} &= l_{T} + b_{T} (\varphi + \varphi^2 +...+ \varphi^h) + S_{T+h-s}
	\end{aligned}
```

As for the t=1,...,s first values there exists no previous values for ``S_t``, 
they have to be estimated before starting the recursion computation. 
[`maSeason`](@ref) will be used for the first ``2*s`` values of the detrended time series. 
Similar for forecast horizions h > s, the same latest s estimated season components
from time periods (T-s-1) to T have to be used over again.


```@docs
holtWinters(y :: Vector, λ₁ :: Real, λ₂ :: Real, λ₃ :: Real, s::Int)
```

### General behavior of holtLinear, brownLinear and holtWinters

By default ``\varphi=1`` in `holtLinear`, `brownLinear` and `holtWinters`
method implementations, so that no damping takes place unless explicitly set.
If no start values ``(l_0, b_0)`` are given for the level and slope parameters, the recursion can
only start in time period 3, as start values will be calculated using the first two observations
as ``l_2 = y_2`` and ``b_2 = y_2 - y_1``. 
Given the input vector y has length N, the output vector f
is still has the length (N+h), as well as the output matrix D still has N rows. The not available
values for periods 1 and 2 are put as NaN[^1] to avoid using a Union{T, Missing} type in julia,
which would make a type recast necessary. 

[^1]:
	> Be aware that Julia follows the IEEE 754 standard for floating-point values.
	> The operation `NaN == NaN` will result in false, which makes e.g. comparisons 
	> of vectors very tricky.

## Optimization

There is no general solution to the problem of selecting the optimal parameter values
for any of the above introduced models. The `Optim.jl` package provides box constrained
algorithms to find the parameters that minimize the sum of squres error, this is equivalent to
using a mean squared error function as loss function. But because of the complexity involved it
is not guarannteed that a global minimum will be found, but rather a local minimum. 

By using this procedure, the internal consistency, which is in favour of achieving the best fit for a
decomposition, is choosen over external validity, which the preferred criterium when making forecasts. 
But the best approach for extensive forecasting is to use a framework or julia package which generally supports cross-validation
for evaluating any forecast model with a horizon h greater than 1. 

### How to get the optimized parameters

By using multiple dispatch the above introduced functions can basically be used without specifying any smoothing
parameters for them to be estimated automatically. The current implementation only allows for all smoothing parameters to be omitted all together. 
In addition, the damping parameter can also be estimated, which is the default, otherwise
it has to be manually set e.g. ``\varphi = 1.0`` for it to have no impact.

As an example, instead of setting the parameters manually like in `holtWinters(data, 0.9, 0.9, 0.9, seasons)`, all
smoothing parameters must be droped for the automatic optimization of the smoothing and dumping parameters and the function simplifies to `holtWinters(data, seasons)`.


A named tuple is returned and the estimated parameters can be accessed in two ways: 

```julia
	f, D, p = holtWinters(data, seasons, φ=1.0)
```

```julia
	res = holtWinters(data, seasons, φ=1.0)
	
	f = res.forecast
	D = res.model
	p = res.parameters
	
```



### Methods for automatic parameter optimization

```@docs
holtWinters(y :: Vector, s::Int)
```


```@docs
holtLinear(y :: Vector)
```


```@docs
brownLinear(y :: Vector)
```

