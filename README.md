TrendDecomposition.jl
=====================

TrendDecomposition.jl is a Julia package for decomposition of time series into trend and cycle components. More generally it provides 
both (stochastic) trend component estimation and forecasting, though not all methods are suitable for forecasting.

**Documentation**: [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://sdbrinkmann.github.io/TrendDecomposition.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://sdbrinkmann.github.io/TrendDecomposition.jl/stable/

Get Started
-----------------
This package is now featured on the official general Julia package registry. 
Simply use Julia's package manager pkg to add TrendDecomposition to your preferred environment.

```Julia
@(v1.11) pkg> add TrendDecomposition

julia> using TrendDecomposition
```

The developing branch of this package can either be employed  by cloning this repository or by using the Julia package manager.
With the package manager simply use the add command:
```Julia
@(v1.11) pkg> add https://github.com/sdBrinkmann/TrendDecomposition.jl
```
> [!IMPORTANT]
> This package is currently under rapid development and follows Semantic Versioning. Until the 1.0.0 release is reached,
>  the API of this package can change with any minor version update, 
> please  consult the documentation of this package after each update when using this package.


