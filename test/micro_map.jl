using Test
using MicroclimateMapper
using MicroclimateMapper: _canonical_data, _is_special_key, _resolve_init, _DEFAULT_INIT
using Microclimate: example_microclimate_problem
using Rasters
using Rasters: X, Y, Ti
using Rasters.Extents: Extent
using Unitful

@testset "_is_special_key" begin
    for k in (:weather, :landcover, :dem, :surface_albedo, :roughness_height)
        @test _is_special_key(k)
    end
    for k in (:vapour_pressure_deficit, :mean_temperature, :foo)
        @test !_is_special_key(k)
    end
end

@testset "_resolve_init" begin
    # nothing → full default NamedTuple with every field nothing
    nt = _resolve_init(nothing, nothing)
    @test nt === _DEFAULT_INIT
    @test nt.soil_temperature === nothing
    @test nt.snow_depth === nothing

    # user NamedTuple wins, missing keys filled from defaults
    user = (; snow_depth = 30.0u"cm")
    merged = _resolve_init(user, nothing)
    @test merged.snow_depth == 30.0u"cm"
    @test merged.soil_temperature === nothing
    @test merged.soil_moisture === nothing
end

@testset "_canonical_data" begin
    rast = Raster(zeros(2, 2), (X(1:2), Y(1:2)))
    data = (
        weather = nothing,
        dem = rast,
        vapour_pressure_deficit = rast,
        mean_temperature = rast,
    )
    canonical = _canonical_data(data)
    @test keys(canonical) == (:vapour_pressure_deficit, :mean_temperature)
    @test canonical.vapour_pressure_deficit === rast
end

@testset "MicroMapModel construction" begin
    # Inner MicroModel uses Microclimate.jl's example builder, so we don't
    # have to wire up the soil profile / hydraulics by hand.
    inner = example_microclimate_problem().model

    model = MicroMapModel(;
        micro_model     = inner,
        dem_source      = SRTM,
        weather_source  = TerraClimate{Historical},
    )

    @test model.micro_model === inner
    @test model.dem_source === SRTM
    @test model.weather_source === TerraClimate{Historical}
    @test model.landcover_source === nothing
    @test model.surface_albedo_source === nothing
    @test model.roughness_height_source === nothing
    @test model.lapse_rate_model isa EnvironmentalLapseRate
    @test length(model.output_layers) == 9
end

@testset "MicroMapProblem construction" begin
    inner = example_microclimate_problem().model
    model = MicroMapModel(;
        micro_model    = inner,
        dem_source     = SRTM,
        weather_source = TerraClimate{Historical},
    )
    area  = Extent(X = (146.0, 146.1), Y = (-36.0, -35.9))
    years = 2000:2000

    problem = MicroMapProblem(; model, area, years, template = SRTM)
    @test problem.model === model
    @test problem.area === area
    @test problem.years === years
    @test problem.template === SRTM
    @test problem.init === nothing
    @test problem.data === (;)

    # Custom init + data overrides
    problem2 = MicroMapProblem(;
        model, area, years, template = SRTM,
        init = (soil_temperature = nothing, soil_moisture = fill(0.1, 10)),
        data = (; vapour_pressure_deficit = Raster(zeros(2, 2), (X(1:2), Y(1:2)))),
    )
    @test problem2.init.soil_moisture == fill(0.1, 10)
    @test haskey(problem2.data, :vapour_pressure_deficit)
end
