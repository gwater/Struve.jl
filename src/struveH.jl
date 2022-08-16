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
        a *= -inv((v + i + T(3)/2) * (i + T(3)/2)) * t2
    end
    return out
end
struveH_power_series_cutoff(nu, x) = x < 5 || nu > evalpoly(x, (-0.75, 0.41, 0.023))

const GAMMA_THREE_HALF(::Type{Float32}) = 0.886226925452758f0
const GAMMA_THREE_HALF(::Type{Float64}) = 0.886226925452758
#const GAMMA_THREE_HALF(::Type{BigFloat}) = big"0.8862269254527580136490837416705725913987747280611935641069038949264556422955125"

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
