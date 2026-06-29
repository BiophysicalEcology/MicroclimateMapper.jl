# ERA5 bindings — hourly global reanalysis from ARCO-ERA5

weather_calendar(::Type{<:ERA5}) = Daily()
native_timestep(::Type{<:ERA5}) = Hourly()
weather_loader(::Type{<:ERA5}) = ContiguousTimeSeries()

function weather_variables(::Type{<:ERA5})
    (
        WeatherVariable(:reference_temperature, :t2m, u"K"),
        WeatherVariable(:eastward_wind, :u10, u"m/s"),
        WeatherVariable(:northward_wind, :v10, u"m/s"),
        WeatherVariable(:dewpoint_temperature, :d2m, u"K"),
        WeatherVariable(:pressure, :sp, u"Pa"),
        # Total cloud cover already 0-1 fraction; no unit, no transform.
        WeatherVariable(:cloud_cover, :tcc, 1),
        # ERA5 radiation is stored as J/(m²·hour) accumulations.
        # Unitful converts that to W/m² on assignment to the canonical buffer.
        WeatherVariable(:global_radiation, :ssrd, u"J/m^2/hr"),
        WeatherVariable(:longwave_radiation, :strd, u"J/m^2/hr"),
        # `:tp` is total precipitation in metres of water per hour; multiply
        # by 1000 to get kg/m² (assuming water density of 1000 kg/m³).
        WeatherVariable(:rainfall, :tp, u"kg/m^2", raw -> raw * 1000.0),
    )
end
