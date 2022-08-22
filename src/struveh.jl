#              Struve Functions of first and second kind
#                     struveh(v, x), struvek(v, x)
#
# The Struve function is computed from its power series [1] for small arguments.
# Large arguments are computed from the asymptotic expansions for struvek(v, x) [2].
# Connection formula using bessely(v, x) [3] are used to compute struveh(v, x) from struvek(v, x) and conversely (struvek -> struveh).
# The power series form [1] can only be used for small x and when nu > x otherwise significant cancellation occurs.
# There is a region where the large argument expansion is not valid and where the power series is prone to cancellation.
# In this region we fall back to computing struvek(v, x) from its integral representation [4] using adaptive gaussian integration (QuadGK.jl).
# TODO: Use a different algorithm in this region instead of relying on quadgk that is more efficient and doesn't allocate.
#
# [1] http://dlmf.nist.gov/11.2.E1
# [2] http://dlmf.nist.gov/11.6.E1
# [3] http://dlmf.nist.gov/11.2.E5
# [4] http://dlmf.nist.gov/11.5.E2
#

"""
struveh(nu, x)

Struve function of the first kind of order `nu`.

External links: [DLMF](http://dlmf.nist.gov/11.2.E1), [Wikipedia](https://en.wikipedia.org/wiki/Struve_function)
"""
struveh(v::Real, x::Real) = _struveh(float(v), float(x))
struveh(v, x) = _H_integral(v, x)

_struveh(v, x::AbstractFloat) = _H_integral(v, x)
_struveh(v::Float16, x::Float16) = Float16(_struveh(Float32(v), Float32(x)))

function _struveh(v, x::T) where T <: Union{Float32, Float64}
    if struvek_large_arg_cutoff(v, x)
        # use http://dlmf.nist.gov/11.2.E6
        return struvek_large_argument(v, x) + bessely(v, x)
    elseif struveh_power_series_cutoff(v, x)
        return struveh_power_series(v, x)
    else
        return _H_integral(v, x)
    end
end

"""
struvek(nu, x)

Struve function of the second kind of order `nu`.

External links: [DLMF](http://dlmf.nist.gov/11.2.E5), [Wikipedia](https://en.wikipedia.org/wiki/Struve_function)
"""
struvek(v::Real, x::Real) = _struvek(float(v), float(x))
struvek(v, x) = _K_integral(v, x)

_struvek(v, x::AbstractFloat) = _K_integral(v, x)
_struvek(v::Float16, x::Float16) = Float16(_struvek(Float32(v), Float32(x)))

function _struvek(v::Real, x::T) where T <: Union{Float32, Float64}
    if struvek_large_arg_cutoff(v, x)
        return struvek_large_argument(v, x)
    elseif struveh_power_series_cutoff(v, x)
        return struveh_power_series(v, x) - bessely(v, x)
    else
        return _K_integral(v, x)
    end
end

# H_{nu}(x) Struve function power series
# http://dlmf.nist.gov/11.2.E1
# can use for x < 5 || nu > -0.75 + 0.41x + 0.023x^2
# struve H can be computed with forward recurrence only when x < nu
# backward recurrence may be stable always?
function struveh_power_series(v, x::T) where T
    MaxIter = 50000
    S = promote_type(T, Float64)
    v, x = S(v), S(x)

    out = zero(S)
    a = (x/2) / (gamma(v + S(3)/2) * gamma(S(3)/2))
    iszero(a) && return a
    t2 = (x/2)^2
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= -inv((v + i + S(3)/2) * (i + S(3)/2)) * t2
    end
    v = (x/2)^(v/2)
    return T((out * v) * v)
end
struveh_power_series_cutoff(nu, x::Float64) = x < 6 || nu > evalpoly(x, (-0.75, 0.41, 0.023))
struveh_power_series_cutoff(nu, x::Float32) = x < 26 || nu > evalpoly(x, (-10.0f0, 0.1f0, 0.012f0))

# K_{nu}(x) using large argument expansion
# http://dlmf.nist.gov/11.6.E1
function struvek_large_argument(v, x::T) where T
    MaxIter = 5000
    S = promote_type(T, Float64)
    v, x = S(v), S(x)

    out = zero(S)
    a = (x/2)^(v-1) * gamma(one(S)/2) / gamma(v + one(S)/2)
    iszero(a) && return a
    t2 = (2/x)^2
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) && break
        a *= (i + one(S)/2)*t2 *  (v - one(S)/2 - i)
    end
    return out / T(pi)
end
struvek_large_arg_cutoff(nu, x) = x > 39 && nu < 1.6*x - 43.0
