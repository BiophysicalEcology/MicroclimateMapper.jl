# gridMET (METDATA) bindings — daily ~4 km gridded weather over the

weather_calendar(::Type{<:GRIDMET}) = Daily()
weather_loader(::Type{<:GRIDMET}) = YearlyTimeSeries()

function weather_variables(::Type{<:GRIDMET})
    (
        WeatherVariable(Temperature(Maximum()), :tmmx, u"K"),
        WeatherVariable(Temperature(Minimum()), :tmmn, u"K"),
        # `:pr` is mm/day depth; numerically equal to kg/m².
        WeatherVariable(Rainfall(), :pr, u"kg/m^2"),
        # gridMET relative humidity is stored as percent — convert to 0-1 fraction.
        WeatherVariable(RelativeHumidity(Minimum()), :rmin, 1, percent_to_fraction),
        WeatherVariable(RelativeHumidity(Maximum()), :rmax, 1, percent_to_fraction),
        WeatherVariable(GlobalRadiation(), :srad, u"W/m^2"),
        # 10 m wind speed — `derive!(:reference_wind_max)` applies the
        WeatherVariable(WindSpeed(), :vs, u"m/s"),
    )
end
