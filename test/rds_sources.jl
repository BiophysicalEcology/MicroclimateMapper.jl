using Test
using Dates
using MicroclimateMapper
using MicroclimateMapper: variables, canonical_name, native_field, fallback_source, fallback_layers,
    weather_calendar, native_timestep, loader, longitude_convention, static_source,
    file_variable, points_load_buffer, Daily, Hourly, MinMax, SixHourly,
    ARCOTimeSeries, MonthlyTimeSeries, YearlyTimeSeries, Longitude180, Longitude360
using RasterDataSources
using RasterDataSources: Land, Static

@testset "ERA5" begin
    for T in (ERA5, ERA5{Google,ARCO}, ERA5{ECMWF,ARCO}, ERA5{ECMWF,ARCO,Land},
              ERA5{ECMWF,OnDemand}, ERA5{ECMWF,OnDemand,Land})
        @test weather_calendar(T) == Daily()
        @test native_timestep(T) == Hourly()
        # Every layer exists in the store
        for var in variables(T)
            @test native_field(var) in RasterDataSources.layers(T)
        end
    end
    @test loader(ERA5) == ARCOTimeSeries()
    @test loader(ERA5{ECMWF,ARCO,Land}) == ARCOTimeSeries()
    @test loader(ERA5{ECMWF,OnDemand}) == MonthlyTimeSeries()
    @test loader(ERA5{ECMWF,OnDemand,Land}) == MonthlyTimeSeries()
    @test longitude_convention(ERA5) == Longitude360()
    @test longitude_convention(ERA5{ECMWF,ARCO}) == Longitude180()
    @test longitude_convention(ERA5{ECMWF,OnDemand,Land}) == Longitude180()
    @test file_variable(ERA5{ECMWF,OnDemand}, :temperature_2m) == :t2m
    @test points_load_buffer(ERA5{ECMWF,OnDemand,Land}) < points_load_buffer(ERA5{ECMWF,OnDemand})

    # ERA5-Land has no cloud cover: it comes from ERA5 in the same store
    @test fallback_source(ERA5{ECMWF,ARCO,Land}) == ERA5{ECMWF,ARCO}
    @test fallback_source(ERA5{ECMWF,OnDemand,Land}) == ERA5{ECMWF,OnDemand}
    @test fallback_layers(ERA5{ECMWF,ARCO,Land}) == (:cloud_cover,)
    @test fallback_source(ERA5{ECMWF,ARCO}) === nothing
end

@testset "BARRA" begin
    T = BARRA{BARRAR2,AUS11}
    @test native_timestep(T) == Hourly()
    @test native_timestep(BARRA{BARRAR2,AUS11,Day}) == MinMax()
    @test static_source(T) == BARRA{BARRAR2,AUS11,Static}
    @test static_source(BARRA{BARRAR2,AUS11,Day}) == BARRA{BARRAR2,AUS11,Static}
    for S in (T, BARRA{BARRAR2,AUS11,Day})
        for var in variables(S)
            layer = native_field(var)
            source = layer === :orog ? static_source(S) : S
            @test layer in RasterDataSources.layers(source)
        end
    end
end

@testset "GRIDMET" begin
    @test static_source(GRIDMET) == GRIDMET{RasterDataSources.Elevation}
    @test file_variable(GRIDMET, :tmmx) == :air_temperature
    @test file_variable(GRIDMET{RasterDataSources.Elevation}, :elev) == :elevation
    for var in variables(GRIDMET)
        layer = native_field(var)
        source = layer === :elev ? static_source(GRIDMET) : GRIDMET
        @test layer in RasterDataSources.layers(source)
    end
end

@testset "NCEP" begin
    T = NCEP{SurfaceFlux,1}
    @test native_timestep(T) == SixHourly()
    @test file_variable(T, :air_2m) == :air
    @test file_variable(T, :uwnd_10m) == :uwnd
    for var in variables(T)
        @test native_field(var) in RasterDataSources.layers(T)
    end
end

@testset "SILO and TerraClimate layers exist" begin
    for T in (SILO, TerraClimate)
        for var in variables(T)
            @test native_field(var) in RasterDataSources.layers(T)
        end
    end
end

@testset "ERA5-Land downloads are de-accumulated" begin
    using Rasters
    using MicroclimateMapper: _from_accumulated
    # Two days at one point, accumulating 1 per hour from 00 UTC
    hours = vcat(fill(24.0, 1), 1.0:23.0, 24.0, 1.0:23.0)
    layer = Raster(reshape(hours, 1, 48), (Dim{:point}([1]), Ti(1:48)))
    hourly = _from_accumulated(ERA5{ECMWF,OnDemand,Land}, :total_precipitation, layer)
    @test all(==(1.0), hourly)
    @test _from_accumulated(ERA5{ECMWF,OnDemand,Land}, :temperature_2m, layer) === layer
    @test _from_accumulated(ERA5{ECMWF,OnDemand}, :total_precipitation, layer) === layer
end
