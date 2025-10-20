# TrendDecomposition.jl

Welcome to the TrendDecomposition.jl documentation.

!!! info
	This is a preliminary version of the documentation.
	The package is also not feature complete until version 1.0,
	thus sometimes there are references to not yet implemented features.

TrendDecomposition.jl is a Julia package for decomposition of time series into trend and cycle components. More generally it provides 
both (stochastic) trend component estimation and forecasting, though not all methods are suitable for forecasting.

By using filters and smoothers the most pragmatic approach to trend decomposition is estimating the trend $t$ and defining
the cyclical component $c$ of time series $y$ as $c = y - t$.
Often it is up to the user of this module to calculate the cyclical components themselves with the computed trend returned from a function 
provided by this module.

The following is a list of already implemented and documented methods:

- Exponential Smoothing
  - Simple exponential smoothing
  - Double exponential smoothing / Brown linear method
  - Holt Linear procedure
  - Holt Winters method
  
- Penalized smoothing
  - Bohlmann Filter / Whittaker-Henderson Smoothing
  - Leser / Hodrick-Prescott (HP) Filter
  - Boosted HP Filter
  
- Moving Average (MA)
- Seasonal Average
- Classical Decomposition by moving averages
