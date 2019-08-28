using Struve
using Test

# using Table 2 from
# https://doi.org/10.1016/S0377-0427(01)00580-5
@testset "H0 zeros" begin
    @test abs(Struve.H0(22.9490276305)) < 1e-10
    @test !(abs(Struve.H0(22.0)) < 1e-10)
end

