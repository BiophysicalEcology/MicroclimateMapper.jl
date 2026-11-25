using Test
using MicroclimateMapper: field_capacity_and_wilting_point
using Unitful

@testset "field capacity and wilting point: hand-computed cross-check" begin
    water_limits = field_capacity_and_wilting_point(20.0, 30.0)
    @test ustrip(u"m^3/m^3", water_limits.field_capacity) ≈ 0.335398 atol=1e-6
    @test ustrip(u"m^3/m^3", water_limits.wilting_point) ≈ 0.173124 atol=1e-6
end
