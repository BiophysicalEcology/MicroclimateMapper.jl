"""
Abstract type for lapse rate formulations. Dispatch on concrete subtypes to select
the temperature adjustment used in mesoclimate corrections.
"""
abstract type LapseRate end

"""
    EnvironmentalLapseRate <: LapseRate

Standard environmental (average) lapse rate: 6.5 K per 1000 m.
"""
struct EnvironmentalLapseRate <: LapseRate end

"""
    DryAdiabaticLapseRate <: LapseRate

Dry adiabatic lapse rate: 9.8 K per 1000 m. Used for unsaturated rising air.
"""
struct DryAdiabaticLapseRate <: LapseRate end

"""
    SaturatedAdiabaticLapseRate <: LapseRate

Saturated (moist) adiabatic lapse rate: ~6.0 K per 1000 m.
Approximates average conditions; varies with temperature and pressure.
"""
struct SaturatedAdiabaticLapseRate <: LapseRate end

"""
    CustomLapseRate{R} <: LapseRate

User-supplied lapse rate. Provide `rate` as a Unitful K/m quantity.

# Example
```julia
lr = CustomLapseRate(0.007u"K/m")
```
"""
struct CustomLapseRate{R} <: LapseRate
    rate::R
end

"""
    lapse_rate(::LapseRate) → Quantity{K/m}

Return the lapse rate constant for the given lapse rate type.
"""
lapse_rate(::EnvironmentalLapseRate)      = 0.0065u"K/m"
lapse_rate(::DryAdiabaticLapseRate)       = 0.0098u"K/m"
lapse_rate(::SaturatedAdiabaticLapseRate) = 0.006u"K/m"
lapse_rate(lr::CustomLapseRate)           = lr.rate

"""
    lapse_adjust_temperature(lr::LapseRate, T, Δz)

Adjust temperature `T` for an elevation difference `Δz` using lapse rate `lr`.

`Δz` = target elevation − reference elevation. Positive `Δz` (target is higher)
lowers temperature; negative `Δz` (target is lower) raises it.

# Examples
```julia
lapse_adjust_temperature(EnvironmentalLapseRate(), 20.0u"°C", 1000.0u"m")  # ≈ 13.5 °C
lapse_adjust_temperature(DryAdiabaticLapseRate(),  20.0u"°C", 1000.0u"m")  # ≈ 10.2 °C
```
"""
lapse_adjust_temperature(lr::LapseRate, T, Δz) = T - lapse_rate(lr) * Δz
