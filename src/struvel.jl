#              Modified Struve Functions of first and second kind
#                     struvel(v, x), struvem(v, x)
#
# The modified struve function is computed from its power series [1] for small arguments.
# Large arguments are computed from the asymptotic expansions for struvem(v, x) [2].
# Connection formula using besseli(v, x) [3] are used to compute struvem(v, x) from struvel(v, x) and conversely (struvel -> struvem).
# The power series form [1] can be computed without cancellation over the whole domain, though it is faster converging for small arguments.
#
# [1] http://dlmf.nist.gov/11.2.E2
# [2] http://dlmf.nist.gov/11.6.E2
# [3] http://dlmf.nist.gov/11.2.E6
#

"""
    struvel(nu, x)

Modified Struve function of the first kind of order `nu`.

External links: [DLMF](http://dlmf.nist.gov/11.2.E2), [Wikipedia](https://en.wikipedia.org/wiki/Struve_function)
"""
struvel(v::Real, x::Real) = _struvel(float(v), float(x))
struvel(v, x) = _L_integral(v, x)

_struvel(v, x::AbstractFloat) = _L_integral(v, x)
_struvel(v::Float16, x::Float16) = Float16(_struvel(Float32(v), Float32(x)))

function _struvel(v, x::T) where T <: Union{Float32, Float64}
    if struvem_large_arg_cutoff(v, x)
        # use http://dlmf.nist.gov/11.2.E6
        return struvem_large_argument(v, x) + besseli(v, x)
    else
        return struvel_power_series(v, x)
    end
end

"""
    struvem(nu, x)

Modified Struve function of the second kind of order `nu`.

External links: [DLMF](http://dlmf.nist.gov/11.2.E6), [Wikipedia](https://en.wikipedia.org/wiki/Struve_function)
"""
struvem(v::Real, x::Real) = _struvem(float(v), float(x))
struvem(v, x) = _M_integral(v, x)

_struvem(v, x::AbstractFloat) = _M_integral(v, x)
_struvem(v::Float16, x::Float16) = Float16(_struvem(Float32(v), Float32(x)))

function _struvem(v, x::T) where T <: Union{Float32, Float64}
    if struvem_large_arg_cutoff(v, x)
        return struvem_large_argument(v, x)
    else
        # use http://dlmf.nist.gov/11.2.E6
        return struvel_power_series(v, x) - besseli(v, x)
    end
end

# L_{nu}(x) Modified Struve function power series
# http://dlmf.nist.gov/11.2.E2
# power series can always be applied without significant cancellation
# could be premature overflow for large orders
# TODO: use either large expansion or use log power series
function struvel_power_series(v, x::T) where T
    MaxIter = 10000
    out = zero(T)
    a = (x/2)^(v+1) / (gamma(v + T(3)/2) * gamma(T(3)/2))
    iszero(a) && return a
    t2 = (x/2)^2
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= inv((v + i + 3/2) * (i + 3/2)) * t2
    end
    return out
end

# M_{nu}(x) using large argument expansion
# http://dlmf.nist.gov/11.6.E2
# valid when nu < 7 - 0.59x + 0.011x^2
function struvem_large_argument(v, x::T) where T
    MaxIter = 50000
    out = zero(T)
    one_half = one(T) / 2
    a = (x/2)^(v-1) * gamma(one_half) / gamma(v + one_half)
    iszero(a) && return a
    t2 = (x/2)^(-2)
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= -(i + one_half) * t2 * (v - one_half - i)
    end
    return out / T(pi)
end

struvem_large_arg_cutoff(nu, x::Float64)  = x > 35 && nu < evalpoly(x, (7.0, -0.59, 0.011))
struvem_large_arg_cutoff(nu, x::Float32)  = x > 32 && nu < 1.5f0*x