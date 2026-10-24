using Test
using Dates
using MicroclimateMapper
using MicroclimateMapper: variables, layers, weather_calendar, native_timestep, loader,
    Daily, MinMax, SingleForecast, _seas5_request_extent, _check_weather_dates,
    _daily_mean, _daily_accumulation, _seas5_read, _load_layers, _SEAS5_FIELDS
using NCDatasets
using RasterDataSources
using Rasters
using Rasters: X, Y, Ti
using Rasters.Extents: Extent

# A file laid out as CDS returns SEAS5 (confirmed against a live download):
# dims (longitude, latitude, forecast_period, forecast_reference_time, number).
function write_cds_seas5(path, f; lons = [144.5], lats = [-37.5], nlead, members = 0:50)
    data = [f(x, y, l, m) for x in lons, y in lats, l in 1:nlead, _ in 1:1, m in members]
    NCDataset(path, "c") do ds
        defDim(ds, "longitude", length(lons)); defDim(ds, "latitude", length(lats))
        defDim(ds, "forecast_period", nlead); defDim(ds, "forecast_reference_time", 1)
        defDim(ds, "number", length(members))
        defVar(ds, "longitude", lons, ("longitude",))
        defVar(ds, "latitude", lats, ("latitude",))
        defVar(ds, "number", collect(members), ("number",))
        defVar(ds, "valid_time", zeros(nlead), ("forecast_period",))
        defVar(ds, "var", data, ("longitude", "latitude", "forecast_period", "forecast_reference_time", "number"))
    end
    return path
end

@testset "SEAS5Forecast" begin
    S = seas5_forecast(Date(2026, 9, 15); member = 3)
    @test S === SEAS5Forecast{Date(2026, 9, 1), 3}     # snapped to the start month
    @test seas5_forecast(Date(2026, 9, 1)) === SEAS5Forecast{Date(2026, 9, 1), 0}
    @test_throws ArgumentError seas5_forecast(Date(2026, 9, 1); member = 51)
    @test weather_calendar(S) == Daily()
    @test native_timestep(S) == MinMax()
    @test loader(S) == SingleForecast()
    @test layers(S) == (:temperature_max, :temperature_min, :dewpoint_temperature,
        :wind_speed, :global_radiation, :rainfall)
    # Every declared native field has a CDS recipe.
    @test all(v -> haskey(_SEAS5_FIELDS, MicroclimateMapper.native_field(v)), variables(S))
end

@testset "SEAS5 request extent" begin
    # Point-like and small areas get the whole-degree cell(s) around them.
    @test _seas5_request_extent(Extent(X = (144.96, 144.96), Y = (-37.81, -37.81))) ==
        Extent(X = (144.0, 145.0), Y = (-38.0, -37.0))
    @test _seas5_request_extent(Extent(X = (144.46, 145.46), Y = (-38.31, -37.31))) ==
        Extent(X = (144.0, 146.0), Y = (-39.0, -37.0))
    @test _seas5_request_extent(Extent(X = (145.0, 145.0), Y = (-37.0, -37.0))) ==
        Extent(X = (145.0, 146.0), Y = (-37.0, -36.0))
    @test _seas5_request_extent(Extent(X = (180.0, 180.0), Y = (90.0, 90.0))) ==
        Extent(X = (179.0, 180.0), Y = (89.0, 90.0))
end

@testset "SEAS5 run dates" begin
    S = seas5_forecast(Date(2026, 9, 1))
    @test _check_weather_dates(S, collect(Date(2026, 9, 1):Day(1):Date(2027, 4, 3))) === nothing  # day 215
    @test_throws ArgumentError _check_weather_dates(S, collect(Date(2026, 8, 31):Day(1):Date(2026, 9, 5)))
    @test_throws ArgumentError _check_weather_dates(S, collect(Date(2027, 4, 1):Day(1):Date(2027, 4, 4)))
    @test _check_weather_dates(TerraClimate{Historical}, [Date(1900, 1, 1)]) === nothing  # default no-op
end

@testset "SEAS5 daily reductions" begin
    six_hourly = reshape(Float64.(1:8), 1, 1, 8)            # two days of 6-hourly values
    @test vec(_daily_mean(six_hourly)) == [2.5, 6.5]
    accumulated = reshape([1.0, 3.0, 3.0 - 1e-7, 6.0], 1, 1, 4)
    @test vec(_daily_accumulation(accumulated)) ≈ [1.0, 2.0, 0.0, 3.0 + 1e-7]  # noise clamped to 0
end

@testset "SEAS5 _seas5_read" begin
    mktempdir() do dir
        f(x, y, l, m) = 1000x + 10l + m
        path = write_cds_seas5(joinpath(dir, "a.nc"), f; lons = [144.5, 145.5], lats = [-37.5, -38.5], nlead = 3)
        r = _seas5_read(path, 7, 3)
        @test r.lons == [144.5, 145.5] && r.lats == [-37.5, -38.5]
        @test r.values == [f(x, y, l, 7) for x in r.lons, y in r.lats, l in 1:3]
        @test_throws ErrorException _seas5_read(path, 7, 4)        # lead count mismatch
        hindcast = write_cds_seas5(joinpath(dir, "h.nc"), f; nlead = 3, members = 0:24)
        @test_throws ErrorException _seas5_read(hindcast, 30, 3)   # member not in file
    end
end

# The whole loader through `getraster`'s cache: files planted where
# `getraster(SEAS5, ...)` would download them, so no CDS key or network.
@testset "SEAS5 _load_layers" begin
    mktempdir() do dir
        withenv("RASTERDATASOURCES_PATH" => dir, "CDSAPI_KEY" => nothing, "HOME" => dir, "USERPROFILE" => dir) do
            init = Date(2026, 9, 1)
            area = Extent(X = (144.96, 144.96), Y = (-37.81, -37.81))
            extent = _seas5_request_extent(area)
            # Per CDS variable: value at lead l (in that file's own steps), member m.
            synthetic = Dict(
                "maximum_2m_temperature_in_the_last_24_hours" => (l, m) -> 290.0 + m,
                "minimum_2m_temperature_in_the_last_24_hours" => (l, m) -> 280.0 + m,
                "2m_dewpoint_temperature" => (l, m) -> 270.0 + mod(l, 4),   # 271,272,273,270 per day
                "10m_u_component_of_wind" => (l, m) -> isodd(l) ? 3.0 : -3.0,
                "10m_v_component_of_wind" => (l, m) -> 4.0,
                "surface_solar_radiation_downwards" => (l, m) -> 1.5e7 * l,  # accumulated
                "total_precipitation" => (l, m) -> 0.002 * l,                # accumulated, m
            )
            for spec in values(_SEAS5_FIELDS), v in spec.cds
                path = RasterDataSources.rasterpath(RasterDataSources.SEAS5, (v,);
                    date = init, extent, leadtime = collect(spec.leadtime), system = "51")
                mkpath(dirname(path))
                g = synthetic[v]
                write_cds_seas5(path, (x, y, l, m) -> g(l, m); nlead = length(spec.leadtime))
            end

            S = seas5_forecast(init; member = 2)
            fields = (:mx2t24, :mn2t24, :d2m, :ws10, :ssrd, :tp)
            stack = _load_layers(SingleForecast(), S, fields, area, 2026:2027)
            sep1 = dayofyear(init)                        # forecast day 1's index
            last_day = sep1 + 214
            @test all(l -> size(l) == (1, 1, 730), values(stack))  # full-years calendar
            @test isnan(stack.mx2t24[1, 1, sep1 - 1]) && isnan(stack.mx2t24[1, 1, last_day + 1])
            @test stack.mx2t24[1, 1, sep1] == 292.0 && stack.mx2t24[1, 1, last_day] == 292.0  # member 2
            @test stack.mn2t24[1, 1, sep1] == 282.0
            @test stack.d2m[1, 1, sep1] ≈ 271.5                  # mean of the day's four 6-hourly values
            @test stack.ws10[1, 1, sep1] ≈ 5.0                   # mean speed, not |mean vector| (4.0)
            @test stack.ssrd[1, 1, sep1 + 10] ≈ 1.5e7            # de-accumulated, J/m² per day
            @test stack.tp[1, 1, sep1 + 10] ≈ 0.002              # m per day (×1000 → 2 kg/m²)
            @test lookup(stack.tp, X) == [144.5] && lookup(stack.tp, Y) == [-37.5]
        end
    end
end
