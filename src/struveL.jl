struveL(v::Real, x::Real) = _struveL(float(v), float(x))
struveL(v, x) = _L_integral(v, x)

_struveL(v, x::AbstractFloat) = _L_integral(v, x)
_struveL(v::Float32, x::Float32) = Float32(_struveL(Float64(v), Float64(x)))
_struveL(v::Float16, x::Float16) = Float16(_struveL(Float32(v), Float32(x)))

function _struveL(v, x::T) where T <: Float64
    if struveM_large_arg_cutoff(v, x)
        # use http://dlmf.nist.gov/11.2.E6
        return struveM_large_argument(v, x) + besseli(v, x)
    else
        return struveL_power_series(v, x)
    end
end

struveM(v::Real, x::Real) = _struveM(float(v), float(x))
struveM(v, x) = _M_integral(v, x)

_struveM(v, x::AbstractFloat) = _M_integral(v, x)
_struveM(v::Float32, x::Float32) = Float32(_struveM(Float64(v), Float64(x)))
_struveM(v::Float16, x::Float16) = Float16(_struveM(Float32(v), Float32(x)))

function _struveM(v, x::T) where T <: Float64
    if struveM_large_arg_cutoff(v, x)
        return struveM_large_argument(v, x)
    else
        # use http://dlmf.nist.gov/11.2.E6
        return struveL_power_series(v, x) - besseli(v, x)
    end
end

# L_{nu}(x) Modified Struve function power series
# http://dlmf.nist.gov/11.2.E2
# power series can always be applied without significant cancellation
# could be premature overflow for large orders
# TODO: use either large expansion or use log power series
function struveL_power_series(v, x::T) where T
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
function struveM_large_argument(v, x::T) where T
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

struveM_large_arg_cutoff(nu, x)  = x > 35 && nu < evalpoly(x, (7.0, -0.59, 0.011))
