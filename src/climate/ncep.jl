# NCEP/NCAR Reanalysis bindings.
#
# Currently supports `NCEP{SurfaceGauss}` — global daily surface-level
# reanalysis on a Gaussian grid, available 1948-present. Each yearly file
# (one per layer per year) contains 365 daily timesteps.
#
# Wind comes as U/V components (`:uwnd_10m`, `:vwnd_10m`), so the chain's
# `derive!(:wind_speed)` step combines them into a scalar speed.
# Humidity comes as specific humidity (`:shum_2m`), so the chain's
# `derive!(:actual_vapour_pressure)` step combines it with surface
# pressure (`:pres`) to recover the vapour pressure the RH derivations
# consume.
#
# `:tcdc` (cloud cover) is only valid for 2005 in the NCEP archive, so we
# don't declare it — `cloud_cover` falls back to the radiation-based
# derivation from `:dswrf`.

temporal_resolution(::Type{<:NCEP}) = DailyResolution()
weather_loader(::Type{<:NCEP{SurfaceGauss}}) = YearlyTimeSeries()

# NCEP getraster requires a `dataset` kwarg ("reanalysis" or "reanalysis2");
# we use the Reanalysis 1 archive throughout.
_extra_getraster_kwargs(::Type{<:NCEP}) = (; dataset = "reanalysis")

function weather_variables(::Type{<:NCEP{SurfaceGauss}})
    (
        WeatherVariable(:maximum_temperature, :tmax, u"K"),
        WeatherVariable(:minimum_temperature, :tmin, u"K"),
        # Wind as horizontal vector components — combined to scalar speed
        # by `derive!(:wind_speed)`.
        WeatherVariable(:u_wind, :uwnd_10m, u"m/s"),
        WeatherVariable(:v_wind, :vwnd_10m, u"m/s"),
        # Specific humidity (kg/kg) and surface pressure (Pa) feed
        # `derive!(:actual_vapour_pressure)`, which then feeds
        # `derive!(:vapour_pressure_deficit)` and the RH derivations.
        WeatherVariable(:specific_humidity, :shum_2m, 1),
        WeatherVariable(:surface_pressure, :pres, u"Pa"),
        WeatherVariable(:downward_shortwave_radiation, :dswrf, u"W/m^2"),
        # `:prate` is precipitation mass flux (kg/m²/s); per day that's
        # `rate × 86400 s` in kg/m². The aggregated daily total then maps
        # straight to `rainfall`.
        WeatherVariable(:rainfall, :prate, u"kg/m^2", raw -> raw * 86400.0),
    )
end
