using Test
using MicroclimateMapper
using MicroclimateMapper: pedotransfer, field_capacity_and_wilting_point,
    CosbyUnivariate, CosbyMultivariate, Campbell1985
using Unitful

@testset "pedotransfer: hand-computed cross-check" begin
    clay, silt, sand, bulk_density = 20.0, 30.0, 50.0, 1.4u"Mg/m^3"

    cosby_uni = pedotransfer(CosbyUnivariate(), clay, silt, sand, bulk_density)
    @test cosby_uni.campbell_b ≈ 5.081399828238082 atol=1e-9
    @test ustrip(u"J/kg", cosby_uni.air_entry_potential) ≈ 1.6458862922770194 atol=1e-9
    @test ustrip(u"kg*s/m^3", cosby_uni.saturated_conductivity) ≈ 0.0005471814266161991 atol=1e-12

    cosby_multi = pedotransfer(CosbyMultivariate(), clay, silt, sand, bulk_density)
    @test cosby_multi.campbell_b ≈ 5.081399828238082 atol=1e-9
    @test ustrip(u"J/kg", cosby_multi.air_entry_potential) ≈ 1.759542771404438 atol=1e-9
    @test ustrip(u"kg*s/m^3", cosby_multi.saturated_conductivity) ≈ 0.0005742901234042081 atol=1e-12

    campbell = pedotransfer(Campbell1985(), clay, silt, sand, bulk_density)
    @test campbell.campbell_b ≈ 6.516628561369164 atol=1e-9
    @test ustrip(u"J/kg", campbell.air_entry_potential) ≈ 2.354775710369497 atol=1e-9
    @test ustrip(u"kg*s/m^3", campbell.saturated_conductivity) ≈ 0.00017701514143601225 atol=1e-12

    fc_pwp = field_capacity_and_wilting_point(clay, silt)
    @test ustrip(u"m^3/m^3", fc_pwp.field_capacity) ≈ 0.335398 atol=1e-6
    @test ustrip(u"m^3/m^3", fc_pwp.wilting_point) ≈ 0.173124 atol=1e-6
end

@testset "pedotransfer: sign convention" begin
    for model in (CosbyUnivariate(), CosbyMultivariate(), Campbell1985())
        for (clay, silt, sand) in ((5.0, 5.0, 90.0), (20.0, 40.0, 40.0), (60.0, 20.0, 20.0))
            result = pedotransfer(model, clay, silt, sand, 1.3u"Mg/m^3")
            @test result.air_entry_potential > 0u"J/kg"
        end
    end
end

@testset "field_capacity_and_wilting_point: physical ordering" begin
    for (clay, silt) in ((5.0, 5.0), (20.0, 40.0), (60.0, 20.0))
        fc_pwp = field_capacity_and_wilting_point(clay, silt)
        @test 0u"m^3/m^3" < fc_pwp.wilting_point < fc_pwp.field_capacity < 1u"m^3/m^3"
    end
end

@testset "pedotransfer: physical sanity bounds across textures" begin
    sandy, loamy, clayey = (5.0, 5.0, 90.0), (20.0, 40.0, 40.0), (60.0, 20.0, 20.0)
    for model in (CosbyUnivariate(), CosbyMultivariate(), Campbell1985())
        for (clay, silt, sand) in (sandy, loamy, clayey)
            result = pedotransfer(model, clay, silt, sand, 1.3u"Mg/m^3")
            @test 2 <= result.campbell_b <= 16  # Campbell1985 runs higher for clay-heavy soils
            @test 1e-8u"kg*s/m^3" < result.saturated_conductivity < 1u"kg*s/m^3"
            @test 0u"J/kg" < result.air_entry_potential < 20u"J/kg"
        end
    end
end
