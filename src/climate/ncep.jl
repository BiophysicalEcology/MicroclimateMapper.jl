# NCEP/NCAR Reanalysis bindings.

weather_calendar(::Type{<:NCEP{<:SurfaceFlux}}) = Daily()
weather_loader(::Type{<:NCEP{<:SurfaceFlux}}) = YearlyTimeSeries()
native_timestep(T::Type{<:NCEP{<:SurfaceFlux}}) = SixHourly()

function weather_variables(::Type{<:NCEP{<:SurfaceFlux}})
    (
        WeatherVariable(Temperature(), :air_2m, u"K"),
        WeatherVariable(EastwardWindSpeed(), :uwnd_10m, u"m/s"),
        WeatherVariable(NorthwardWindSpeed(), :vwnd_10m, u"m/s"),
        WeatherVariable(SpecificHumidity(), :shum_2m, 1),
        WeatherVariable(Pressure(), :pres, u"Pa"),
        WeatherVariable(GlobalRadiation(), :dswrf, u"W/m^2"),
        WeatherVariable(LongwaveRadiation(), :dlwrf, u"W/m^2"),
        # prate is kg/m²/s; × 6 h × 3600 s/h → kg/m² per 6 h block.
        WeatherVariable(Rainfall(), :prate, u"kg/m^2", raw -> raw * 6 * 3600),
    )
end
