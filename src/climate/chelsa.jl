# CHELSA bindings: `CHELSA{Climate}` (1981-2010 monthly climatology) and
# `CHELSA{Future{Climate, ...}}`

# ---------------------------------------------------------------------------
# CHELSA{Climate} — historical baseline climatology
# ---------------------------------------------------------------------------

weather_loader(::Type{<:CHELSA{Climate}}) = MonthlyClimatology()

function weather_variables(::Type{<:CHELSA{Climate}})
    (
        WeatherVariable(Temperature(Maximum()), :tasmax, u"°C"),
        WeatherVariable(Temperature(Minimum()), :tasmin, u"°C"),
        WeatherVariable(WindSpeed(), :sfcWind, u"m/s"),
        WeatherVariable(GlobalRadiation(), :rsds, u"W/m^2"),
        WeatherVariable(Rainfall(), :pr, u"kg/m^2"),
        # CHELSA gives a single mean relative humidity — use it for both
        # reference_humidity min and max so the vapour-pressure-deficit-based
        # derivation is skipped on both ends.
        WeatherVariable(RelativeHumidity(Minimum()), :hurs, 1, percent_to_fraction),
        WeatherVariable(RelativeHumidity(Maximum()), :hurs, 1, percent_to_fraction),
        # `cloud_cover` is native, so the radiation-based cloud derivation
        # is skipped; cloud_min/max derivations still fire from this value.
        WeatherVariable(CloudCover(), :clt, 1, percent_to_fraction),
    )
end

# ---------------------------------------------------------------------------
# CHELSA{Future{Climate, ...}} — projected period climatology
# ---------------------------------------------------------------------------

weather_loader(::Type{<:CHELSA{<:Future{Climate}}}) = MonthlyClimatologyPeriod()

primary_layers(::Type{<:CHELSA{<:Future{Climate}}}) =
    (:temperature_max, :temperature_min, :rainfall)
fallback_layers(::Type{<:CHELSA{<:Future{Climate}}}) =
    (:wind_speed, :global_radiation, :reference_humidity_min, :reference_humidity_max, :cloud_cover)
fallback_source(::Type{<:CHELSA{<:Future{Climate}}}) = CHELSA{Climate}

function weather_variables(::Type{<:CHELSA{<:Future{Climate}}})
    (
        WeatherVariable(Temperature(Maximum()), :tmax, u"°C"),
        WeatherVariable(Temperature(Minimum()), :tmin, u"°C"),
        WeatherVariable(Rainfall(), :prec, u"kg/m^2"),
        # Fallback fields from the CHELSA{Climate} 1981-2010 baseline:
        WeatherVariable(WindSpeed(), :sfcWind, u"m/s"),
        WeatherVariable(GlobalRadiation(), :rsds, u"W/m^2"),
        WeatherVariable(RelativeHumidity(Minimum()), :hurs, 1, percent_to_fraction),
        WeatherVariable(RelativeHumidity(Maximum()), :hurs, 1, percent_to_fraction),
        WeatherVariable(CloudCover(), :clt, 1, percent_to_fraction),
    )
end
