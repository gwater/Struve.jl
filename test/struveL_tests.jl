#
# Tests for struveL and struve M
#

# implement test based on http://dlmf.nist.gov/11.4.E7
for x in 0.1:1.1:300.0
    @test isapprox(Struve.struveL.(0.5, x), sqrt(2 / (BigFloat(x)*pi)) * (cosh(BigFloat(x)) - 1), rtol=1e-14)
end

# implement test based on http://dlmf.nist.gov/11.4.E11
for x in 0.1:1.1:300.0
    s = Struve.struveL.(3/2, x)
    x = BigFloat(x)
    t = -sqrt(x / pi / 2) * (1 - 2/x^2) + sqrt(2 / x / pi) * (sinh(x) - cosh(x) / x)
    @test isapprox(s, t, rtol=1e-14)
end

for x in rand(100)*200, v in rand(100)*100
    t = Struve.struveL_power_series(BigFloat(v), BigFloat(x))
    @test isapprox(Struve.struveL(v, x), t, rtol=1e-13)
end

# testing for struveM which can be prone to cancellation due to subtraction
# in some domains struveL ≈ besseli leading to loss of digits
for x in 0.1:0.15:2.0, v in 0.1:1.2:20.0
    x = x*v
    t = Struve.struveL_power_series(BigFloat(v), BigFloat(x)) - besseli(v, x)
    @test isapprox(Struve.struveM(v, x), t, rtol=1e-9)
end
