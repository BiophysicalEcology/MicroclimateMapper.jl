using Test

@testset "MicroclimateMapper" begin
    @testset "landcover"      begin include("landcover.jl") end
    @testset "terrain"        begin include("terrain.jl") end
    @testset "output writing" begin include("output_writing.jl") end
    @testset "micro raster"   begin include("micro_raster.jl") end
    @testset "micro vector"   begin include("micro_vector.jl") end
    @testset "soil water limits" begin include("soil_water_limits.jl") end
end
