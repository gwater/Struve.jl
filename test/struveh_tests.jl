#
# Tests for struveh and struvek
#

# implement test based on http://dlmf.nist.gov/11.4.E5
for x in 0.1:1.5:100.0
    @test isapprox(Struve.struveh.(0.5, x), sqrt(2 / (BigFloat(x)*pi)) * (1 - cos(BigFloat(x))), rtol=1e-11)
end

# implement test based on http://dlmf.nist.gov/11.4.E9
for x in 0.1:1.1:150.0
    s = Struve.struveh.(3/2, x)
    x = BigFloat(x)
    t = sqrt(x / pi / 2) * (1 + 2/x^2) - sqrt(2 / x / pi) * (sin(x) + cos(x) / x)
    @test isapprox(s, t, rtol=1e-12)
end

for x in 0.1:1.3:80.0, v in 0.1:1.2:40.0
    t = Struve.struveh_power_series(BigFloat(v), BigFloat(x))
    @test isapprox(Struve.struveh(v, x), t, rtol=1e-10)
end

# need to increase precision for BigFloat calculations
setprecision(BigFloat, 2500)
# test orders between 100-500
for x in [0.3, 0.5, 0.8, 0.9, 1.0, 1.1, 1.5, 2.0, 3.0], v in [100, 120.2, 170.0, 200.32, 250.2, 300.2, 400.5, 450.0, 490.5]
    x = x*v
    t = Struve.struveh_power_series(BigFloat(v), BigFloat(x))
    @test isapprox(Struve.struveh(v, x), t, rtol=1e-10)
end
# test orders between 1000-2000
for x in [0.9, 1.0, 1.05, 1.1, 1.2], v in [1000.0, 1200.2, 1500.2, 1752.2345, 2000.0]
    x = x*v
    t = Struve.struveh_power_series(BigFloat(v), BigFloat(x))
    @test isapprox(Struve.struveh(v, x), Float64(t), rtol=1e-10)
end

# testing for struvek
for x in 0.1:0.15:2.0, v in 0.1:1.2:40.0
    x = x*v
    t = Struve.struveh_power_series(BigFloat(v), BigFloat(x)) - bessely(v, x)
    @test isapprox(Struve.struvek(v, x), t, rtol=1e-9)
end

for x in [0.75, 0.9, 1.0, 1.1, 1.5, 2.0, 3.0], v in [100.0, 223.2, 400.0, 600.2, 1000.0, 1200.25, 2000.0]
    x = x*v
    t = Struve.struveh_power_series(BigFloat(v), BigFloat(x)) - bessely(v, x)
    @test isapprox(Struve.struvek(v, x), Float64(t), rtol=1e-9)
end
