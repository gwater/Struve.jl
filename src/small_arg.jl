# H_{nu}(x) Struve function power series
# http://dlmf.nist.gov/11.2.E1
# can use for x < 5 || nu > -0.75 + 0.41x + 0.023x^2
# could result in premature underflow for nu >> x
# struve H can be computed with forward recurrence only when x < nu
# backward recurrence may be stable always?
function struveH_power_series(v, x::T) where T
    MaxIter = 50000
    out = zero(T)
    a = (x/2)^(v+1) / (gamma(v + T(3)/2) * GAMMA_THREE_HALF(T))
    iszero(a) && return a
    t2 = (x/2)^2
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= -inv((v + i + 3/2) * (i + 3/2)) * t2
    end
    return out
end
struveH_power_series_cutoff(nu, x) = x < 5 || nu > evalpoly(x, (-0.75, 0.41, 0.023))

# L_{nu}(x) Modified Struve function power series
# http://dlmf.nist.gov/11.2.E2
# power series can always be applied without significant cancellation
function struveL_power_series(v, x::T) where T
    MaxIter = 50000
    out = zero(T)
    a = (x/2)^(v+1) / (gamma(v + T(3)/2) * GAMMA_THREE_HALF(T))
    iszero(a) && return a
    t2 = (x/2)^2
    for i in 0:MaxIter
        out += a
        abs(a) < eps(T) * abs(out) && break
        a *= inv((v + i + 3/2) * (i + 3/2)) * t2
    end
    return out
end

const GAMMA_THREE_HALF(::Type{Float32}) = 0.886226925452758f0
const GAMMA_THREE_HALF(::Type{Float64}) = 0.886226925452758
#const GAMMA_THREE_HALF(::Type{BigFloat}) = big"0.8862269254527580136490837416705725913987747280611935641069038949264556422955125"

function struveK(nu, x::T) where T
    if struveH_power_series_cutoff(nu, x)
        return struveH_power_series(nu, x) - bessely(nu, x)
    else
        return struveK_large_argument(nu, x)
    end
end

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

function struveM(nu, x::T) where T
    if x > 5*nu && x > 20.0
        return struveM_large_argument(nu, x)
    else
        return struveL_power_series(nu, x) - besseli(nu, x)
    end
end

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
