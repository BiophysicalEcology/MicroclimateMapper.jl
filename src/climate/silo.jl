# SILO (Scientific Information for Land Owners) bindings — daily ~5 km
# gridded weather over Australia, 1889-present. One annual NetCDF per
# variable per year, downloaded and cached via `RasterDataSources.SILO`.
#
# SILO has no native wind field — falls back to CRUCL2's monthly wind
# climatology (`:wnd`), same as NicheMapR's own SILO/AWAP implementations.
#
# `:rh_tmax` (RH at the time of Tmax, the day's lowest) maps to
# `reference_humidity_min`; `:rh_tmin` (RH at Tmin, the day's highest) maps
# to `reference_humidity_max` — same convention as GRIDMET's `:rmin`/`:rmax`.

temporal_resolution(::Type{<:SILO}) = DailyResolution()
weather_loader(::Type{<:SILO}) = YearlyTimeSeries()

primary_layers(::Type{<:SILO}) = (:max_temp, :min_temp, :daily_rain, :rh_tmax, :rh_tmin, :radiation)
fallback_layers(::Type{<:SILO}) = (:wnd,)
fallback_source(::Type{<:SILO}) = CRUCL2

function weather_variables(::Type{<:SILO})
    (
        WeatherVariable(:maximum_temperature, :max_temp, u"°C"),
        WeatherVariable(:minimum_temperature, :min_temp, u"°C"),
        # `:daily_rain` is mm/day depth; numerically equal to kg/m².
        WeatherVariable(:rainfall, :daily_rain, u"kg/m^2"),
        WeatherVariable(:reference_humidity_min, :rh_tmax, 1, _silo_percent_to_fraction),
        WeatherVariable(:reference_humidity_max, :rh_tmin, 1, _silo_percent_to_fraction),
        # Daily total solar exposure (MJ/(m²·day)); Unitful converts to
        # W/m² on assignment to the canonical buffer.
        WeatherVariable(:downward_shortwave_radiation, :radiation, u"MJ/m^2/d"),
        # Fallback: CRUCL2's monthly wind climatology, expanded to daily.
        WeatherVariable(:wind_speed, :wnd, u"m/s"),
    )
end

_silo_percent_to_fraction(raw) = raw / 100.0
