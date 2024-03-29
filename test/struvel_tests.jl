#
# Tests for struvel and struve M
#

# implement test based on http://dlmf.nist.gov/11.4.E7
for x in 0.1:1.1:300.0
    @test isapprox(Struve.struvel.(0.5, x), sqrt(2 / (BigFloat(x)*pi)) * (cosh(BigFloat(x)) - 1), rtol=1e-14)
end

# implement test based on http://dlmf.nist.gov/11.4.E11
for x in 0.1:1.1:300.0
    s = Struve.struvel.(3/2, x)
    x = BigFloat(x)
    t = -sqrt(x / pi / 2) * (1 - 2/x^2) + sqrt(2 / x / pi) * (sinh(x) - cosh(x) / x)
    @test isapprox(s, t, rtol=1e-14)
end

for x in 0.1:1.5:100.0, v in 0.0:0.3:60.0
    t = Float64(Struve.struvel_power_series(BigFloat(v), BigFloat(x)))
    @test isapprox(Struve.struvel(v, x), t, rtol=1e-12)
end

# testing for struvem which can be prone to cancellation due to subtraction
# in some domains struvel ≈ besseli leading to loss of digits
for x in 0.1:0.15:2.0, v in 0.1:1.2:20.0
    x = x*v
    t = Struve.struvel_power_series(BigFloat(v), BigFloat(x)) - besseli(v, x)
    @test isapprox(Struve.struvem(v, x), t, rtol=1e-9)
end
