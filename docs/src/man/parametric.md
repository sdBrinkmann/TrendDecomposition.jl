# Parametric Estimation

In general, parametric models in regards to time series analysis refer to (mixed) autoregressive and moving
average process models. Here, they are denoted in following way:

```math
	y_{t} = \sum_{i=1}^{p} \phi_{i} y_{t-i} + \sum_{i=1}^{q} \theta_{i} \epsilon_{t-i}  + \epsilon_{t},
```

where ``\phi_{i}`` are the parameters of an autoregressive process of p-th order and ``\theta_{i}``
are the parameters of an moving average process of q-th order and ``\epsilon`` denotes the error term. 

Most of the following functions are internally used by other estimators in this package, but
because of their usefullness, they can be directly called  by the user. Since for now wrapper functions
are missing for each estimation method there exists a seperate function. For the special cases q = 0, use the `ar` functions and for p = 0 use `maNR`.

## Autoregressiv Process AR(p)

```@docs
arBurg
```

```@docs
arOLS
```

```@docs
arYuleWalker
```

```@docs
arDurbinLevinson
```

## Moving Average Process MA(p)
A moving average process, abbreviated as MA(q) can be estimated by maximum likelihood using the Newton-Raphson. This implementation follows Anderson[^1], for that reason  the parameters ``\theta`` are denoted here as ``\alpha``.
```@docs
maNR
```
## ARMA(p, q)

For now mixed autoregressive and moving average process, abbreviated as ARMA(p, q) can be estimated by maximum likelihood using the Newton-Raphson
procedure. Both p and q have to be greater than 0.
This estimation method was proposed by Åström and Torsten[^2] and specifically this implementation follows 
the method proposed by Anderson[^1], for that reason the parameter ``\phi`` are denoted here as ``\beta`` and the parameters ``\theta`` as ``\alpha``,
as in the original paper. Normally, the algorithm converges very fast within 4 to 8 iterations, but sometimes the algorithm becomes unstable which could lead in this implementation to a singular matrix execption. In that case, choosing different starting points for the parameters can help to find a stable solution.

```@docs
armaNR
```


[^1]:
	> Anderson, T. W. (1977). Estimation for autoregressive moving average models in the time and frequency domains. The Annals of Statistics, 5(5), 842-865.
	
[^2]:
	> Åström, K. J., & Torsten, B. (1965). Numerical identification of linear dynamic systems from normal operating records. IFAC Proceedings Volumes, 2(2), 96-111.
