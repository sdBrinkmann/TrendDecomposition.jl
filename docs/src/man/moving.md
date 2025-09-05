# Moving Average (MA)

## Rolling Average
The name rolling average instead of moving average is choosen here, because by
default the functions works like a rolling window that slides through the 
entire time series from ``t = 1,...,T``, calculating a value for each datum and 
it is up to the user to deceide to have the boundary values discarded. 

For centered MA the choosen order p has to be odd and thus can be written
as ``p = 2k + 1`` so that

```math
	MA_t(p) = \frac{1}{p} \sum_{i=-k}^{k} x_{t+i}
```
With `centered=false` the last (p-1) lagging values are used instead, for this option the
order p can be either even or odd. 

```math
	MA_t(p) = \frac{1}{p} \sum_{i=0}^{p-1} x_{t-i}
```

```@docs
rollingAverage(y :: Vector, p :: Int; centered::Bool = true, discard::Bool = false, offset::Int = 0)
```

Both methods allow an offset for calculating the MA. With a positive offset an additional lead value is included at
the expense of a lagged value and the reverse holds for negative offsets. Negativ offset are only allowed for 
centered MA, since for `centered=false` the maximum amount of lagged values are already used by default.

### MA with weights

One can also specify weights ``\omega_i`` for each i-th term, subject to that they sum to unity, e.g. for centered MA

```math
	WMA_t(p) = \frac{1}{p} \sum_{i=-k}^{k} \omega_i x_{t+i} \qquad  s.t \sum_{i=-k}^{k} \omega_{i} = 1.
```

For a MA of order p, the weights provided by the user have to be a vector of size (p x 1).

```@docs
rollingAverage(y :: Vector, p :: Int, weights :: Vector;centered::Bool = true, discard::Bool = false, offset::Int = 0)
```

### Calculation for boundary or edge values

When using a centered moving average of order p, for the first ``\frac{p-1}{2}`` and the last ``\frac{p-1}{2}`` values 
of the time series there are not enought data points available. To generalise about these boundary cases, the normal
centered MA can also be thought as a special case of weighted MA with uniform weights, where each weight must equal 
``w_i = \frac{1}{p}``. For k missing values, the weights that still can be applied are normalized so that they
sum to 1, so that each of the (p-k) remaining weights now equals ``\frac{1}{p-k}``. 
As an example, the calculation of the first value of a time series y using a centered MA of order 5 becomes 
``MA_{1}(5) = \sum_{i=1}^{3} \frac{1}{3} y_i = \frac{1}{3} \sum_{i=1}^{3} y_i``.

For weighted MA with given weights ``w_1,..,w_p``, in case that for the first k weights there are no data points available for the 
calculation, the remaining weights are normalized as follows ``w'_i = \frac{w_i}{\sum_{j=k+1}^{p}w_j}``. Correspondingly,
in case the last k weights cannot be used, the remaining ``p-k`` weights are normalized: ``w'_i = \frac{w_i}{\sum_{j=1}^{p-k}w_j}``

As it is common to drop data for which the full order of the specified MA cannot be computed, 
`discard=true` will return NaN[^1] for each datum in the output vector, where
not enough data points are available due to boundary cases.

[^1]:
	> Be aware that Julia follows the IEEE 754 standard for floating-point values.
	> The operation `NaN == NaN` will result in false, which makes e.g. comparisons 
	> of vectors very tricky.

## Seasonal Average
For each season component ``S_i`` this function calculates the average for all observation that where recorded in the same type of season given ``1,2,...s`` types.

Since the formal mathematical notation requires the use of sets, the interested reader
can see the formula for calculating the seasonal averages in the notes down below[^2].

[^2]:
	> For the time series ``y_1, y_2,...,y_T`` of length T and the number of seasons s, each observation ``y_t, t \in \{1,...,T\}``, 
	> can belong only to one set ``\mathcal{S_i}``, ``i \in \{1,..,s\}`` out of all possible s sets, where each set is comprised of all observations that take place in the same type of season. 
	> Therefore the superset ``\mathcal{S}``, consisting of all individual disjoint sets ``\mathcal{S_i}``, can be defined as  
	> ```math 
	> \mathcal{S} = \bigcup_{i = 1}^{s} \mathcal{S}_i,
	> ```
	> where  ``\mathcal{S}_i = \{y_t | ((t-1) \mod s) + 1 = i\}`` and ``\mathcal{S_i} \cap \mathcal{S_j} = \emptyset`` for all ``i \neq j``.
	> So the seasonal average for the i-th number of season ``S_i`` can be written and computed as
	> ```math 
	> S_i = \frac{1}{|\mathcal{S}_i|} \sum_{x \in \mathcal{S}_i} x
	>
	> ```
	


```@docs
maSeason
```

## Classical decomposition

Given the implementations of [`maSeason`](@ref) and [`rollingAverage`](@ref),
the following two decompositions using moving averages can be archived: 

The additive model
```math 
	y_t = \tau_t + S_t + u_t
```
or the multiplicative model
```math 
	y_t = \tau_t * S_t * u_t,
```
where ``\tau_t`` is the trend component estimated using moving averages,
``S_t`` is the average season component estimated using the detrended time series and
``u_t`` is simply the residual factor defined as ``u_t = y_t - \tau_t - S_t``.

```@docs
maDecompose
```
