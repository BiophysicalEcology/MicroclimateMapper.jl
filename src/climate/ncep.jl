# NCEP/NCAR Reanalysis bindings.

weather_calendar(::Type{<:NCEP{<:SurfaceFlux}}) = Daily()
weather_loader(::Type{<:NCEP{<:SurfaceFlux}}) = YearlyTimeSeries()
native_timestep(T::Type{<:NCEP{<:SurfaceFlux}}) = SixHourly()

function weather_variables(::Type{<:NCEP{<:SurfaceFlux}})
    (
        WeatherVariable(:reference_temperature, :air_2m, u"K"),
        WeatherVariable(:eastward_wind, :uwnd_10m, u"m/s"),
        WeatherVariable(:northward_wind, :vwnd_10m, u"m/s"),
        WeatherVariable(:specific_humidity, :shum_2m, 1),
        WeatherVariable(:pressure, :pres, u"Pa"),
        WeatherVariable(:global_radiation, :dswrf, u"W/m^2"),
        WeatherVariable(:longwave_radiation, :dlwrf, u"W/m^2"),
        # prate is kg/m²/s; × 6 h × 3600 s/h → kg/m² per 6 h block.
        WeatherVariable(:rainfall, :prate, u"kg/m^2", raw -> raw * 6 * 3600),
    )
end
