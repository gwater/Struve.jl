
# based on https://dlmf.nist.gov/11.13.iii, works for real and complex arguments
# fails at z=0 because Struve K and Bessel Y diverge there

integrate(f, a, b) = QuadGK.quadgk(f, a, b)[1]

# only valid for real(z) > 0, https://dlmf.nist.gov/11.5.E2

_K0_integral(z::T) where T = (2 / T(pi)) * integrate(t -> exp(-z * sinh(t)), zero(T), T(Inf))

function _K_integral(ν, z::T) where T
    out = 2 * (z / 2)^ν / (sqrt(T(π)) * gamma(ν + one(T)/2))
    return out * integrate(t -> exp(-z * t) * (one(T) + t^2)^(ν - one(T)/2), zero(T), T(Inf))
end

_H0_integral(z::Real) = _K0_integral(z) + bessely0(z)
_H_integral(ν, z::Real) = _K_integral(ν, z) + bessely(ν, z)

# only valid for real(ν) > -0.5, https://dlmf.nist.gov/11.5.E1
function _H0_integral(z::Complex{T}) where {T <: Real}
    return (2 / T(π)) * integrate(ϑ -> sin(z * cos(ϑ)), zero(T), T(π)/2)
end

function _H_integral(ν, z::Complex{T}) where {T <: Real} 
    out = 2 * (z / 2)^ν / (sqrt(T(π)) * gamma(ν + one(T)/2))
    return out * integrate(t -> (1 - t^2)^(ν - one(T)/2) * sin(z * t), zero(T), one(T))
end

_K0_integral(z::Complex{T}) where {T <: Real} = _H0_integral(z) - bessely0(z)
_K_integral(ν, z::Complex{T}) where {T <: Real} = _H_integral(ν, z) - bessely(ν, z)

# only valid for real(ν) > -0.5, https://dlmf.nist.gov/11.5.E4
_M0_integral(z) = -(2 / T(π)) * integrate(ϑ -> exp(-z * cos(ϑ)), zero(T), T(π)/2)


function _M_integral(ν, z)
    out =  -2 * (z / 2)^ν / (sqrt(T(π)) * gamma(ν + one(T)/2))
    return out * integrate(t -> exp(-z * t) * (1 - t^2)^(ν - one(T)/2), zero(T), one(T))
end

_L0_integral(z)   = besseli(0, z) + _M0_integral(z)
_L_integral(ν, z) = besseli(ν, z) + _M_integral(ν, z)
