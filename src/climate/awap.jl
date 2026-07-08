# AWAP (Australian Water Availability Project) bindings.
#
# Daily ~5 km gridded Australian weather, 1900-present. One file per layer
# per day (`getraster(AWAP, layer; date)` returns a single `.grid` path).
#
# AWAP provides actual vapour pressure at 09:00 (`:vprpress09`) and 15:00
# (`:vprpress15`). For now we use only the 09:00 value as the day's
# `actual_vapour_pressure` — that's a simplification (the daily mean would
# need a multi-input merge step we don't yet have). The afternoon value
# is left unused.
#
# AWAP has no native wind field — falls back to CRUCL2's monthly wind
# climatology (`:wnd`), same as NicheMapR's own SILO/AWAP implementations.

temporal_resolution(::Type{<:AWAP}) = DailyResolution()
weather_loader(::Type{<:AWAP}) = DailyFiles()

primary_layers(::Type{<:AWAP}) = (:tmax, :tmin, :rainfall, :solar, :vprpress09)
fallback_layers(::Type{<:AWAP}) = (:wnd,)
fallback_source(::Type{<:AWAP}) = CRUCL2

function weather_variables(::Type{<:AWAP})
    (
        WeatherVariable(:maximum_temperature, :tmax, u"°C"),
        WeatherVariable(:minimum_temperature, :tmin, u"°C"),
        # `:rainfall` is mm/day depth; numerically equal to kg/m^2.
        WeatherVariable(:rainfall, :rainfall, u"kg/m^2"),
        # AWAP solar is MJ/(m^2·day); Unitful converts to W/m^2 on
        # assignment to the canonical buffer.
        WeatherVariable(:downward_shortwave_radiation, :solar, u"MJ/m^2/d"),
        # 09:00 vapour pressure in hPa — converted to kPa on assignment.
        # TODO: when a multi-input merge step exists, average with
        # `:vprpress15` for a true daily-mean actual vapour pressure.
        WeatherVariable(:actual_vapour_pressure, :vprpress09, u"hPa"),
        # Fallback: CRUCL2's monthly wind climatology, expanded to daily.
        WeatherVariable(:wind_speed, :wnd, u"m/s"),
    )
end
