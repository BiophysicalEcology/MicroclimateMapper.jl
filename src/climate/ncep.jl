# NCEP/NCAR Reanalysis bindings.
#
# `NCEP` is RasterDataSources' data-source type `NCEP{Group, Reanalysis, Period}`
# — a *kind of data*, with time on its own `Period` axis. MicroclimateMapper
# never wraps it in a time-flavoured type; it reads the kind from the data
# source and the timestep separately.
#
# We bind the `SurfaceFlux` group (the T62 Gaussian grid, ~1.875°): 6-hourly
# 2 m temperature, 10 m wind, 2 m specific humidity, surface pressure,
# short/longwave radiation, and precipitation in one consistent group. This
# avoids the much larger pressure-level files (`PressureLevels`) which carry
# all 17 standard levels when only the surface-adjacent level is needed.
#
# Time:
#   * calendar = Daily (every day; the solver runs consecutive-day mode).
#   * native timestep = translated from the source's `Period` parameter via
#     `RasterDataSources.period` — six-hourly when the parameter is omitted.
#     The generic `resample!` pass then carries each quantity to the run's
#     target timestep (hourly by default): linear interpolation for the met
#     fields, solar-geometry disaggregation for shortwave, block-sum for
#     rainfall.
#
# RDS serves SurfaceFlux natively at six-hourly (1460 steps/year) for every
# layer; `_load_weather` errors if a layer arrives at anything else.

weather_calendar(::Type{<:NCEP{<:SurfaceFlux}}) = Daily()
weather_loader(::Type{<:NCEP{<:SurfaceFlux}}) = YearlyTimeSeries()

# Translation layer: the source's native `Period` (six-hourly if the third
# type parameter is omitted) → MicroclimateMapper's within-day timestep.
native_timestep(T::Type{<:NCEP{<:SurfaceFlux}}) = SixHourly()

# Wind arrives as U/V components and humidity as specific humidity + surface
# pressure; the assembler's shared combiners turn these into scalar wind speed
# and vapour pressure at the native timestep before resampling. Names are the
# canonical native-tier quantities.
function weather_variables(::Type{<:NCEP{<:SurfaceFlux}})
    (
        WeatherVariable(:reference_temperature,        :air_2m,   u"K"),
        WeatherVariable(:u_wind,                       :uwnd_10m, u"m/s"),
        WeatherVariable(:v_wind,                       :vwnd_10m, u"m/s"),
        WeatherVariable(:specific_humidity,            :shum_2m,  1),
        WeatherVariable(:pressure,                     :pres,     u"Pa"),
        WeatherVariable(:global_radiation,             :dswrf,    u"W/m^2"),
        WeatherVariable(:longwave_radiation,           :dlwrf,    u"W/m^2"),
        # prate is kg/m²/s; × 6 h × 3600 s/h → kg/m² per 6 h block.
        WeatherVariable(:rainfall,                     :prate,    u"kg/m^2",
                        raw -> raw * 6 * 3600),
    )
end
