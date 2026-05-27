using Test

@testset "MicroclimateMapper" begin
    @testset "landcover"      begin include("landcover.jl") end
    @testset "terrain"        begin include("terrain.jl") end
    @testset "output writing" begin include("output_writing.jl") end
    @testset "micro map"      begin include("micro_map.jl") end
    @testset "micro points"   begin include("micro_points.jl") end
end
