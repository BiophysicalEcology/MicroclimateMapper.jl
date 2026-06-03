using Test
using Dates
using MicroclimateMapper
using MicroclimateMapper: _make_points_dim, _points_extent,
    _to_points, _to_points_constant, _to_points_with_time,
    _stack_to_points
using Microclimate: example_microclimate_problem, example_soil_profile
using Rasters
using Rasters: X, Y, Ti, Near
using Rasters.DimensionalData: Dim, MergedLookup, hasdim
using Unitful

const POINTS = [(146.05, -35.95), (146.02, -35.92), (146.08, -35.98)]

@testset "_make_points_dim" begin
    pdim = _make_points_dim(POINTS)
    @test pdim isa Dim{:point}
    @test length(pdim) == length(POINTS)
    lkp = lookup(pdim)
    @test lkp isa MergedLookup
    @test collect(lkp) == [(146.05, -35.95), (146.02, -35.92), (146.08, -35.98)]
end

@testset "_points_extent" begin
    ext = _points_extent(POINTS)
    @test ext.X == (146.02, 146.08)
    @test ext.Y == (-35.98, -35.92)
end

@testset "_to_points scalar" begin
    pdim = _make_points_dim(POINTS)
    r = _to_points(0.15, pdim)
    @test size(r) == (3,)
    @test all(==(0.15), r)
    @test dims(r) == (pdim,)
end

@testset "_to_points 2-D raster" begin
    pdim = _make_points_dim(POINTS)
    # Synthetic 2-D raster: value(x, y) = x + 100 * y
    xs = 146.00:0.01:146.10
    ys = -36.00:0.01:-35.90
    data = [x + 100 * y for x in xs, y in ys]
    rast = Raster(data, (X(xs), Y(ys)))

    out = _to_points(rast, pdim)
    @test size(out) == (3,)
    @test dims(out) == (pdim,)
    # Each point gets the nearest cell's value
    for (i, (x, y)) in enumerate(POINTS)
        @test out[i] == rast[X(Near(x)), Y(Near(y))]
    end
end

@testset "_to_points 3-D raster preserves Ti" begin
    pdim = _make_points_dim(POINTS)
    xs = 146.00:0.01:146.10
    ys = -36.00:0.01:-35.90
    nt = 4
    data = [x + 100 * y + 0.01 * t for x in xs, y in ys, t in 1:nt]
    rast = Raster(data, (X(xs), Y(ys), Ti(1:nt)))

    out = _to_points(rast, pdim)
    @test size(out) == (3, nt)
    @test hasdim(out, Ti)
    @test hasdim(out, Dim{:point})
    for (i, (x, y)) in enumerate(POINTS), t in 1:nt
        @test out[i, t] == rast[X(Near(x)), Y(Near(y)), Ti(t)]
    end
end

@testset "_to_points preserves units" begin
    pdim = _make_points_dim(POINTS)
    xs = 146.00:0.01:146.10
    ys = -36.00:0.01:-35.90
    data = [(x + 100 * y) * u"°C" for x in xs, y in ys]
    rast = Raster(data, (X(xs), Y(ys)))
    out = _to_points(rast, pdim)
    @test eltype(out) === eltype(rast)
end

@testset "_stack_to_points" begin
    pdim = _make_points_dim(POINTS)
    xs = 146.00:0.01:146.10
    ys = -36.00:0.01:-35.90
    a = Raster([x for x in xs, _ in ys], (X(xs), Y(ys)); name = :a)
    b = Raster([y for _ in xs, y in ys], (X(xs), Y(ys)); name = :b)
    stack = RasterStack((; a, b))
    out = _stack_to_points(stack, pdim)
    @test keys(out) == (:a, :b)
    @test dims(out.a) == (pdim,)
    @test dims(out.b) == (pdim,)
end

@testset "MicroVectorProblem construction" begin
    inner = example_microclimate_problem().model
    model = MicroMapModel(;
        micro_model = inner,
        dem_source = SRTM,
        weather_source = TerraClimate{Historical},
    )
    dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)
    soil_profile = example_soil_profile(inner.depths)

    problem = MicroVectorProblem(; model, points = POINTS, dates, soil_profile)
    @test problem.model === model
    @test problem.points === POINTS
    @test problem.dates === dates
    @test problem.soil_profile === soil_profile
    @test problem.init === nothing
    @test problem.data === (;)

    # Single-day run
    problem_day = MicroVectorProblem(;
        model, points = POINTS, dates = Date(2000, 6, 29), soil_profile)
    @test problem_day.dates === Date(2000, 6, 29)

    problem2 = MicroVectorProblem(;
        model, points = POINTS, dates, soil_profile,
        init = (soil_temperature = nothing, soil_moisture = fill(0.1, 10)),
        data = (; vapour_pressure_deficit = Raster(zeros(2, 2), (X(1:2), Y(1:2)))),
    )
    @test problem2.init.soil_moisture == fill(0.1, 10)
    @test haskey(problem2.data, :vapour_pressure_deficit)
end
