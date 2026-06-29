# WorldClim bindings

# ---------------------------------------------------------------------------
# WorldClim{Climate} — 1970-2000 baseline climatology
# ---------------------------------------------------------------------------

weather_loader(::Type{<:WorldClim{Climate}}) = MonthlyClimatology()

function weather_variables(::Type{<:WorldClim{Climate}})
    (
        WeatherVariable(Temperature(Maximum()), :tmax, u"°C"),
        WeatherVariable(Temperature(Minimum()), :tmin, u"°C"),
        WeatherVariable(WindSpeed(), :wind, u"m/s"),
        # WorldClim radiation is kJ/(m^2·day); Unitful converts to W/m^2 on
        # assignment to the canonical buffer (`Quantity{W/m^2}`).
        WeatherVariable(GlobalRadiation(), :srad, u"kJ/m^2/d"),
        # `:prec` is mm depth; 1 mm of rainfall ≡ 1 kg/m^2.
        WeatherVariable(Rainfall(), :prec, u"kg/m^2"),
        # `:vapr` is actual vapour pressure — feeds the
        # vapour-pressure-deficit derivation that runs before relative humidity.
        WeatherVariable(ActualVapourPressure(), :vapr, u"kPa"),
    )
end

# ---------------------------------------------------------------------------
# WorldClim{Future{Climate, ...}} — projected period climatology, 20-year windows.
# ---------------------------------------------------------------------------

weather_loader(::Type{<:WorldClim{<:Future{Climate}}}) = MultiBandClimatologyPeriod()

primary_layers(::Type{<:WorldClim{<:Future{Climate}}}) =
    (:temperature_max, :temperature_min, :rainfall)
fallback_layers(::Type{<:WorldClim{<:Future{Climate}}}) =
    (:wind_speed, :global_radiation, :actual_vapour_pressure)
fallback_source(::Type{<:WorldClim{<:Future{Climate}}}) = WorldClim{Climate}

function weather_variables(::Type{<:WorldClim{<:Future{Climate}}})
    (
        WeatherVariable(Temperature(Maximum()), :tmax, u"°C"),
        WeatherVariable(Temperature(Minimum()), :tmin, u"°C"),
        WeatherVariable(Rainfall(), :prec, u"kg/m^2"),
    )
end
