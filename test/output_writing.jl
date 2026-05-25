using Test
using BiophysicalGrids
using BiophysicalGrids: _write_slice!
using Rasters
using Rasters: X, Y, Ti, Dim

@testset "_write_slice! 3D scalar layer" begin
    nx, ny, nt = 4, 3, 5
    rast = Raster(zeros(nx, ny, nt), (X(1:nx), Y(1:ny), Ti(1:nt)))
    src = collect(10.0 .* (1:nt))

    # (X(2), Y(3)) singles out one (X, Y) pixel; all `Ti` values fill in.
    _write_slice!(rast, src, (X(2), Y(3)))
    @test rast[X=2, Y=3, Ti=1] == 10.0
    @test rast[X=2, Y=3, Ti=5] == 50.0
    @test rast[X=1, Y=1, Ti=1] == 0.0  # untouched
end

@testset "_write_slice! 4D layer (soil-style)" begin
    nx, ny, nt, nd = 3, 2, 4, 6
    rast = Raster(zeros(nx, ny, nt, nd),
        (X(1:nx), Y(1:ny), Ti(1:nt), Dim{:depth}(1:nd)))
    src = [10.0*t + d for t in 1:nt, d in 1:nd]

    _write_slice!(rast, src, (X(2), Y(1)))
    @test rast[X(2), Y(1), Ti(3), Dim{:depth}(4)] == 34.0
    @test rast[X(2), Y(1), Ti(1), Dim{:depth}(1)] == 11.0
    @test rast[X(1), Y(1), Ti(1), Dim{:depth}(1)] == 0.0  # untouched
end

@testset "_write_slice! 4D layer (profile-style, :height dim)" begin
    # The :height variant tests that the trailing-dim detection in
    # `_write_slice!` works regardless of whether the dim is `:depth` or
    # `:height` — both kinds of 4D output should write to the right cell.
    nx, ny, nt, nh = 3, 2, 4, 2
    rast = Raster(zeros(nx, ny, nt, nh),
        (X(1:nx), Y(1:ny), Ti(1:nt), Dim{:height}(1:nh)))
    src = [100.0*h + t for t in 1:nt, h in 1:nh]

    _write_slice!(rast, src, (X(3), Y(2)))
    @test rast[X(3), Y(2), Ti(4), Dim{:height}(2)] == 204.0
    @test rast[X(3), Y(2), Ti(1), Dim{:height}(1)] == 101.0
end
