# gridMET (METDATA) bindings — daily ~4 km gridded weather over the

weather_calendar(::Type{<:GRIDMET}) = Daily()
weather_loader(::Type{<:GRIDMET}) = YearlyTimeSeries()

function weather_variables(::Type{<:GRIDMET})
    (
        WeatherVariable(:maximum_temperature, :tmmx, u"K"),
        WeatherVariable(:minimum_temperature, :tmmn, u"K"),
        # `:pr` is mm/day depth; numerically equal to kg/m².
        WeatherVariable(:rainfall, :pr, u"kg/m^2"),
        # gridMET relative humidity is stored as percent — convert to 0-1 fraction.
        WeatherVariable(:reference_humidity_min, :rmin, 1, percent_to_fraction),
        WeatherVariable(:reference_humidity_max, :rmax, 1, percent_to_fraction),
        WeatherVariable(:global_radiation, :srad, u"W/m^2"),
        # 10 m wind speed — `derive!(:reference_wind_max)` applies the
        WeatherVariable(:wind_speed, :vs, u"m/s"),
    )
end
