using Test

@testset "BiophysicalGrids" begin
    @testset "landcover"      begin include("landcover.jl") end
    @testset "terrain"        begin include("terrain.jl") end
    @testset "output writing" begin include("output_writing.jl") end
end
