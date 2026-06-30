# WorldClim bindings

# ---------------------------------------------------------------------------
# WorldClim{Climate} — 1970-2000 baseline climatology
# ---------------------------------------------------------------------------

loader(::Type{<:WorldClim{Climate}}) = MonthlyClimatology()

function variables(::Type{<:WorldClim{Climate}})
    (
        Variable(Temperature(Maximum()), :tmax, u"°C"),
        Variable(Temperature(Minimum()), :tmin, u"°C"),
        Variable(WindSpeed(), :wind, u"m/s"),
        # WorldClim radiation is kJ/(m^2·day); Unitful converts to W/m^2 on
        # assignment to the canonical buffer (`Quantity{W/m^2}`).
        Variable(GlobalRadiation(), :srad, u"kJ/m^2/d"),
        # `:prec` is mm depth; 1 mm of rainfall ≡ 1 kg/m^2.
        Variable(Rainfall(), :prec, u"kg/m^2"),
        # `:vapr` is actual vapour pressure — feeds the
        # vapour-pressure-deficit derivation that runs before relative humidity.
        Variable(ActualVapourPressure(), :vapr, u"kPa"),
    )
end

# ---------------------------------------------------------------------------
# WorldClim{Future{Climate, ...}} — projected period climatology, 20-year windows.
# ---------------------------------------------------------------------------

loader(::Type{<:WorldClim{<:Future{Climate}}}) = MultiBandClimatologyPeriod()
fallback_source(::Type{<:WorldClim{<:Future{Climate}}}) = WorldClim{Climate}

function variables(::Type{<:WorldClim{<:Future{Climate}}})
    (
        Variable(Temperature(Maximum()), :tmax, u"°C"),
        Variable(Temperature(Minimum()), :tmin, u"°C"),
        Variable(Rainfall(), :prec, u"kg/m^2"),
    )
end
