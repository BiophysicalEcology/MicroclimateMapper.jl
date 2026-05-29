# WorldClim bindings: `WorldClim{Climate}` (1970-2000 monthly climatology)
# and `WorldClim{Future{Climate, CMIP6, Model, Scenario}}` (projected period
# climatology, 20-year windows).
#
# WorldClim provides `:vapr` (actual vapour pressure in kPa) instead of
# `:vpd`, so it relies on the default chain's
# `vapour_pressure_deficit ← e_s(mean_T) - actual_vapour_pressure` step
# (which fires automatically since VPD isn't in the variable map).
#
# WorldClim{Future{Climate}} only ships `:tmin`, `:tmax`, `:prec` (packed
# as a single multi-band GeoTIFF, 12 monthly bands per file) so it falls
# back to WorldClim{Climate} for wind, radiation, and vapour pressure.

# ---------------------------------------------------------------------------
# WorldClim{Climate} — 1970-2000 baseline climatology
# ---------------------------------------------------------------------------

weather_loader(::Type{<:WorldClim{Climate}}) = MonthlyClimatology()

function weather_variables(::Type{<:WorldClim{Climate}})
    (
        WeatherVariable(:maximum_temperature,          :tmax, u"°C"),
        WeatherVariable(:minimum_temperature,          :tmin, u"°C"),
        WeatherVariable(:wind_speed,                   :wind, u"m/s"),
        # WorldClim radiation is kJ/(m^2·day); Unitful converts to W/m^2 on
        # assignment to the canonical buffer (`Quantity{W/m^2}`).
        WeatherVariable(:downward_shortwave_radiation, :srad, u"kJ/m^2/d"),
        # `:prec` is mm depth; 1 mm of rainfall ≡ 1 kg/m^2.
        WeatherVariable(:rainfall,                     :prec, u"kg/m^2"),
        # `:vapr` is actual vapour pressure — feeds the VPD-from-actual
        # derivation that runs before RH.
        WeatherVariable(:actual_vapour_pressure,       :vapr, u"kPa"),
    )
end

# ---------------------------------------------------------------------------
# WorldClim{Future{Climate, ...}} — projected period climatology
# ---------------------------------------------------------------------------

weather_loader(::Type{<:WorldClim{<:Future{Climate}}})     = MultiBandFutureClimatology()

primary_layers(::Type{<:WorldClim{<:Future{Climate}}})  = (:tmin, :tmax, :prec)
fallback_layers(::Type{<:WorldClim{<:Future{Climate}}}) = (:wind, :srad, :vapr)
fallback_source(::Type{<:WorldClim{<:Future{Climate}}}) = WorldClim{Climate}

function weather_variables(::Type{<:WorldClim{<:Future{Climate}}})
    (
        WeatherVariable(:maximum_temperature,          :tmax, u"°C"),
        WeatherVariable(:minimum_temperature,          :tmin, u"°C"),
        WeatherVariable(:rainfall,                     :prec, u"kg/m^2"),
        # Fallback fields from the WorldClim{Climate} 1970-2000 baseline:
        WeatherVariable(:wind_speed,                   :wind, u"m/s"),
        WeatherVariable(:downward_shortwave_radiation, :srad, u"kJ/m^2/d"),
        WeatherVariable(:actual_vapour_pressure,       :vapr, u"kPa"),
    )
end
