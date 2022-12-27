using Sad
using Test

@testset "Cross-sections" begin
    r_xs = Sad.Rectangular(100., 5., 2., 0.00001, 0.03)
    d_xs = Sad.Dingman(100., 5., 2., 2.5, 0.00001, 0.03)
    @test width(r_xs) == 100.0
    @test depth(r_xs) == 2.0
    @test width(d_xs) ≈ 69.314 atol=0.001
    @test depth(d_xs) ≈ 1.428 atol=0.001
end
