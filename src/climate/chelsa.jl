# CHELSA bindings: `CHELSA{Climate}` (1981-2010 monthly climatology) and
# `CHELSA{Future{Climate, ...}}`

# ---------------------------------------------------------------------------
# CHELSA{Climate} — historical baseline climatology
# ---------------------------------------------------------------------------

weather_loader(::Type{<:CHELSA{Climate}}) = MonthlyClimatology()

function weather_variables(::Type{<:CHELSA{Climate}})
    (
        WeatherVariable(:maximum_temperature, :tasmax, u"°C"),
        WeatherVariable(:minimum_temperature, :tasmin, u"°C"),
        WeatherVariable(:wind_speed, :sfcWind, u"m/s"),
        WeatherVariable(:global_radiation, :rsds, u"W/m^2"),
        WeatherVariable(:rainfall, :pr, u"kg/m^2"),
        # CHELSA gives a single mean relative humidity — use it for both
        # reference_humidity min and max so the vapour-pressure-deficit-based
        # derivation is skipped on both ends.
        WeatherVariable(:reference_humidity_min, :hurs, 1, percent_to_fraction),
        WeatherVariable(:reference_humidity_max, :hurs, 1, percent_to_fraction),
        # `cloud_cover` is native, so the radiation-based cloud derivation
        # is skipped; cloud_min/max derivations still fire from this value.
        WeatherVariable(:cloud_cover, :clt, 1, percent_to_fraction),
    )
end

# ---------------------------------------------------------------------------
# CHELSA{Future{Climate, ...}} — projected period climatology
# ---------------------------------------------------------------------------

weather_loader(::Type{<:CHELSA{<:Future{Climate}}}) = MonthlyClimatologyPeriod()

primary_layers(::Type{<:CHELSA{<:Future{Climate}}}) = (:tmin, :tmax, :prec)
fallback_layers(::Type{<:CHELSA{<:Future{Climate}}}) = (:sfcWind, :hurs, :clt, :rsds)
fallback_source(::Type{<:CHELSA{<:Future{Climate}}}) = CHELSA{Climate}

function weather_variables(::Type{<:CHELSA{<:Future{Climate}}})
    (
        WeatherVariable(:maximum_temperature, :tmax, u"°C"),
        WeatherVariable(:minimum_temperature, :tmin, u"°C"),
        WeatherVariable(:rainfall, :prec, u"kg/m^2"),
        # Fallback fields from the CHELSA{Climate} 1981-2010 baseline:
        WeatherVariable(:wind_speed, :sfcWind, u"m/s"),
        WeatherVariable(:global_radiation, :rsds, u"W/m^2"),
        WeatherVariable(:reference_humidity_min, :hurs, 1, percent_to_fraction),
        WeatherVariable(:reference_humidity_max, :hurs, 1, percent_to_fraction),
        WeatherVariable(:cloud_cover, :clt, 1, percent_to_fraction),
    )
end
