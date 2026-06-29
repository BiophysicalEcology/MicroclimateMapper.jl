# AWAP (Australian Water Availability Project) bindings.
#
# Daily ~5 km gridded Australian weather, 1900-present. One file per layer
# per day (`getraster(AWAP, layer; date)` returns a single `.grid` path).
#
# AWAP measures actual vapour pressure twice a day, at 09:00 and 15:00. Both are
# declared as timed samples; the within-day vapour-pressure curve and the derived
# relative humidity follow generically from that.

weather_calendar(::Type{<:AWAP}) = Daily()
weather_loader(::Type{<:AWAP}) = DailyFiles()

function weather_variables(::Type{<:AWAP})
    (
        WeatherVariable(Temperature(Maximum()), :tmax, u"°C"),
        WeatherVariable(Temperature(Minimum()), :tmin, u"°C"),
        # `:rainfall` is mm/day depth; numerically equal to kg/m^2.
        WeatherVariable(Rainfall(), :rainfall, u"kg/m^2"),
        # AWAP solar is MJ/(m^2·day); Unitful converts to W/m^2 on assignment.
        WeatherVariable(GlobalRadiation(), :solar, u"MJ/m^2/d"),
        # Morning and afternoon vapour-pressure readings in hPa (→ kPa on assignment).
        WeatherVariable(ActualVapourPressure(ClockTime(9)), :vprpress09, u"hPa"),
        WeatherVariable(ActualVapourPressure(ClockTime(15)), :vprpress15, u"hPa"),
    )
end
