
struveH(v::Real, x::Real) = _struveH(float(v), float(x))
struveH(v, x) = _H_integral(v, x)

_struveH(v, x::AbstractFloat) = _H_integral(v, x)
_struveH(v::Float32, x::Float32) = Float32(_struveH(Float64(v), Float64(x)))
_struveH(v::Float16, x::Float16) = Float16(_struveH(Float32(v), Float32(x)))

function _struveH(v, x::T) where T <: Float64
    if struveK_large_arg_cutoff(v, x)
        # use http://dlmf.nist.gov/11.2.E6
        return struveK_large_argument(v, x) + bessely(v, x)
    elseif struveH_power_series_cutoff(v, x)
        return struveH_power_series(v, x)
    else
        return _H_integral(v, x)
    end
end

struveK(v::Real, x::Real) = _struveK(float(v), float(x))
struveK(v, x) = _K_integral(v, x)

_struveK(v, x::AbstractFloat) = _K_integral(v, x)
_struveK(v::Float32, x::Float32) = Float32(_struveK(Float64(v), Float64(x)))
_struveK(v::Float16, x::Float16) = Float16(_struveK(Float32(v), Float32(x)))

function _struveK(v::Real, x::T) where T <: Float64
    if struveK_large_arg_cutoff(v, x)
        return struvek_large_argument(v, x)
    elseif struveH_power_series_cutoff(v, x)
        return struveH_power_series(v, x) - bessely(v, x)
    else
        return _K_integral(v, x)
    end
end

# H_{nu}(x) Struve function power series
# http://dlmf.nist.gov/11.2.E1
# can use for x < 5 || nu > -0.75 + 0.41x + 0.023x^2
# could result in premature underflow for nu >> x
# struve H can be computed with forward recurrence only when x < nu
# backward recurrence may be stable always?
function struveH_power_series(v, x::T) where T
    MaxIter = 50000
    out = zero(T)
    a = (x/2)^(v+1) / (gamma(v + T(3)/2) * gamma(T(3)/2))
    iszero(a) && return a
    t2 = (x/2)^2
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= -inv((v + i + T(3)/2) * (i + T(3)/2)) * t2
    end
    return out
end
struveH_power_series_cutoff(nu, x) = x < 6 || nu > evalpoly(x, (-0.75, 0.41, 0.023))

# K_{nu}(x) using large argument expansion
# http://dlmf.nist.gov/11.6.E1
function struveK_large_argument(v, x::T) where T
    MaxIter = 500
    out = zero(T)
    a = (x/2)^(v-1) * gamma(T(1)/2) / gamma(v + 1/2)
    iszero(a) && return a
    t2 = (x/2)^(-2)
    for i in 0:MaxIter
        out += a
        abs(a) < 1e-14 && break
        a *= (i + 1/2)*t2 *  (v - 1/2 - i)
    end
    return out / T(pi)
end
struveK_large_arg_cutoff(nu, x) = x > 35 && nu < 1.6*x - 39.0
