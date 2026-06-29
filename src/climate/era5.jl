# ERA5 bindings — hourly global reanalysis from ARCO-ERA5

weather_calendar(::Type{<:ERA5}) = Daily()
native_timestep(::Type{<:ERA5}) = Hourly()
weather_loader(::Type{<:ERA5}) = ContiguousTimeSeries()

function weather_variables(::Type{<:ERA5})
    (
        WeatherVariable(Temperature(), :t2m, u"K"),
        WeatherVariable(EastwardWindSpeed(), :u10, u"m/s"),
        WeatherVariable(NorthwardWindSpeed(), :v10, u"m/s"),
        WeatherVariable(DewpointTemperature(), :d2m, u"K"),
        WeatherVariable(Pressure(), :sp, u"Pa"),
        # Total cloud cover already 0-1 fraction; no unit, no transform.
        WeatherVariable(CloudCover(), :tcc, 1),
        WeatherVariable(GlobalRadiation(), :ssrd, u"J/m^2/hr"),
        WeatherVariable(LongwaveRadiation(), :strd, u"J/m^2/hr"),
        # `:tp` is total precipitation in metres of water per hour; multiply
        # by 1000 to get kg/m² (assuming water density of 1000 kg/m³).
        WeatherVariable(Rainfall(), :tp, u"kg/m^2", raw -> raw * 1000.0),
    )
end
