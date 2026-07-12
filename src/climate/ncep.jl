# NCEP/NCAR Reanalysis bindings.

weather_calendar(::Type{<:NCEP{<:SurfaceFlux}}) = Daily()
loader(::Type{<:NCEP{<:SurfaceFlux}}) = YearlyTimeSeries()
native_timestep(T::Type{<:NCEP{<:SurfaceFlux}}) = SixHourly()

function variables(::Type{<:NCEP{<:SurfaceFlux}})
    (
        Variable(Reference(Temperature()), :air_2m, u"K"),
        Variable(EastwardWindSpeed(), :uwnd_10m, u"m/s"),
        Variable(NorthwardWindSpeed(), :vwnd_10m, u"m/s"),
        Variable(SpecificHumidity(), :shum_2m, 1),
        Variable(Pressure(), :pres, u"Pa"),
        Variable(GlobalRadiation(), :dswrf, u"W/m^2"),
        Variable(LongwaveRadiation(), :dlwrf, u"W/m^2"),
        # prate is kg/m²/s; × 6 h × 3600 s/h → kg/m² per 6 h block.
        Variable(Rainfall(), :prate, u"kg/m^2", raw -> raw * 6 * 3600),
    )
end
