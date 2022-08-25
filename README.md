# Struve.jl

This package provides methods to compute the
[Struve functions](https://dlmf.nist.gov/11) H, K, L, and M.

[![Build Status](https://travis-ci.org/gwater/Struve.jl.svg?branch=master)](https://travis-ci.org/gwater/Struve.jl)
[![codecov](https://codecov.io/gh/gwater/Struve.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gwater/Struve.jl)

<img alt="plot of four Struve functions (H₀, K₀, L₀, M₀) on the real axis"
src="./example.png" width="480">

The default methods currently use a mixture of power series, large argument expansions, and integral representations of the functions
which are computed numerically using [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) and [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl).



`Struve.jl` exports four functions:

```julia
struveh(ν, z)
struvek(ν, z)
sruvel(ν, z)
struvem(ν, z)
```
which compute the Struve functions of the first and second kind (`struveh` and `struvek`), and the modified Struve functions of the first and second kind (`struvel` and `struvem`).

It also implements
[fast approximations for H₀ and H₁](http://dx.doi.org/10.1121/1.4968792) on the
real axis (with absolute error below 2×10⁻³).
For fast, high accuracy approximations [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl) may be used.
```julia
Struve.H0_fast(x)
Struve.H1_fast(x)
```

Please note: Implementations are verified against function's power series computed in higher precision as well as explicit forms. In general, we try to provide maximum relative errors throughout the entire range better than `1e-13` in double precision. There may be some regions around cutoffs that are less accurate. Bug reports are very welcome as it is difficult to test over all ranges of order and argument. Rigorous implementations of Struve functions are difficult to find but comparisons to Mathematica or [Keisan](https://keisan.casio.com/exec/system/1222676451) would be helpful. The power series can also be called directly (e.g., `Struve.struveh_power_series(big"1.5", big"90.0")`) using higher precision which will be a good comparison for the lower precisions calculations.
