# Spectral Analysis


## Spectral Density AR(p) process
	
Given that an autoregressive process ``X`` of p'th order is defined in a
regressive fashion as 
```math
X_t = a_1 X_{t-1} + a_2 X_{t-2} + ... + a_p X_{t-p} + \epsilon_t
```
the two functions `arSpectrum` compute an estimate of the spectral density 
of the above defined AR(p) process. 

Either the already estimated or predetermined parameters with coefficents (``a_1, a_2, ..., a_p``) and
variance of the error terms ``\sigma_{\epsilon}`` can be used or for a given time series ``y_t`` the parameters can
be estimated in a first step with one of the selected estimators (Burg, OLS, Yule-Walker).

For both methods the spectral density function is calculated as a function ``h``
of the frequency component ``\omega`` as 
```math
 h(\omega) = \frac{\sigma_{\epsilon}}{2\pi} \frac{1}{|1 - a_1 exp(-i\omega)
 - ... - a_p exp(-ip\omega)|^2}
```

For an early treatment of this procedure see Akaike (1969)[^1]

[^1]:
	> Akaike, H. (1969). Power spectrum estimation through autoregressive model fitting. Annals of the institute of Statistical Mathematics, 21(1), 407-419.
	
	

```@docs
arSpectrum
```

## Periodogram
For a time series ``X(t)`` the periodogram is defined as the modulus-squared of the finite Fourier transform in accordance with the definition of Schulz (1898)[^2]:

[^2]:
	>Schuster, A. (1898), On the investigation of hidden periodicities with application to a supposed 26 day period of meteorological phenomena, Terr. Magn., 3(1), 13–41, doi:10.1029/TM003i001p00013.

```math
	I(\omega) = \frac{1}{2 \pi T} | \sum_{t=0}^{T-1} exp(-i\omega t) X(t)|^2,
```

for ``-\pi \leq \omega \leq \pi`` and ``\omega_j = \frac{2 \pi j}{T}``, (``j = 1,...,m`` where ``m = \lceil (T-1) / 2 \rceil``) 

With regard to the autocovariance function ``\gamma (p)`` one can write the periodogram also 
as 

```math
  I(\omega) = \frac{1}{ 2\pi}  \sum_{j=-T-1}^{T-1} \gamma(j) cos(j \omega).
```

Using the factor ``\frac{1}{2 \pi}`` here lets the periodogram also serve as a direct estimate of
the spectral density or power spectrum of a covariance stationary process (see Brillinger 1975[^3]). 

[^3]:
	>Brillinger, D. R. (1975). Time series: data analysis and theory. Holt, Rinehart and Winston, New York


```@docs
periodogram
```

## Band-Pass Filter

	
	
### Baxter-King 
Baxter-King (1999)[^4] provide an approximate band-pass filter, denoted by the authors as ``BK_{K}(pl, pu)``, where K states the numbers of lags (or leads) to include and the bandwidth is stated in terms of the periodicity of cycles (``2\pi / \omega``) to include (pl denotes the lowest and pu the highest period included).

[^4]:
	> Baxter, M., & King, R. G. (1999). Measuring business cycles: approximate band-pass filters for economic time series. Review of economics and statistics, 81(4), 575-593.
```@docs
baxterKing
```
