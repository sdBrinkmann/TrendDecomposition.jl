# Stochastic trend analysis

## Berveridge-Nelson Decomposition

The Beveridge-Nelson decomposition[^1] divides a time series into two stochastic
components, given that the time series is non-stationary (has a unit root) or
in other terms that the series can be modeled as an ARIMA(p, 1, q) process with drift. 

[^1]:
	> Beveridge, S., & Nelson, C. R. (1981). A new approach to decomposition of economic time series into permanent and transitory components with particular attention to measurement of the ‘business cycle’. Journal of Monetary economics, 7(2), 151-174.

Unlike in other models the decomposition of time series $y$ in terms of trend $\tau$ and cycle $c$ components is defined as:

```math
	y_t = \tau_t - c_t
```



This package offers two ways of computing the decomposition, one by Miller[^2]
for autoregressive models and one by Newbold[^3] for ARMA models. For the latter
it is also possible to calculate the decomposition with already estimated parameters of
an ARIMA(p, 1, q) model, otherwise the function `armaNR` (see [ARMA(p, q)](@ref)) is used.

[^2]:
	> Miller, S. M. (1988). The Beveridge-Nelson decomposition of economic time series: Another economical computational method. Journal of Monetary Economics, 21(1), 141-142.

[^3]:
	> Newbold, P. (1990). Precise and efficient computation of the Beveridge-Nelson decomposition of economic time series. Journal of Monetary Economics, 26(3), 453-457.


```@docs
beveridgeNelson
```
