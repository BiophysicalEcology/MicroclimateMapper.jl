using Test
using Dates
using MicroclimateMapper
using MicroclimateMapper: _check_weather_covers, _ti_range_for_dates, weather_calendar, native_timestep
using RasterDataSources
using Rasters
using Rasters: X, Y, Ti

# A stack with `n` steps on its Ti axis, plus a static layer (which has none).
stack_with(n) = RasterStack((;
    t2m = Raster(zeros(1, 1, n), (X([144.5]), Y([-37.5]), Ti(1:n))),
    elevation = Raster(zeros(1, 1), (X([144.5]), Y([-37.5]))),
))

@testset "_check_weather_covers" begin
    # ERA5-Land's store ending on 2026-08-29 (241 days of hourly data), with
    # a run that reaches 2026-08-31: the case that used to BoundsError.
    years = 2026:2026
    dates = collect(Date(2026, 6, 1):Day(1):Date(2026, 8, 31))
    _, ti_end, _ = _ti_range_for_dates(weather_calendar(ERA5ECMWFLand), native_timestep(ERA5ECMWFLand), years, dates)
    short = stack_with(241 * 24)
    err = try
        _check_weather_covers(ERA5ECMWFLand, short, years, ti_end)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("only available to 2026-08-29", err.msg)

    # Data reaching the last run date (or beyond) passes.
    @test _check_weather_covers(ERA5ECMWFLand, stack_with(ti_end), years, ti_end) === nothing
    @test _check_weather_covers(ERA5ECMWFLand, stack_with(365 * 24), years, ti_end) === nothing

    # Monthly sources count months, not days.
    @test_throws "only 8 of the 12 time steps" _check_weather_covers(TerraClimate, stack_with(8), years, 12)
end
