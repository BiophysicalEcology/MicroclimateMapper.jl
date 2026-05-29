# TerraClimate bindings: a yearly NetCDF time series (one file per year,
# each containing 12 months). Files supply temperature min/max, wind,
# vapour pressure deficit, shortwave radiation, precipitation, and soil
# moisture. The scenario (Historical, Plus2C, Plus4C) is encoded in the
# `TerraClimate{Scenario}` source-type parameter.

# TerraClimate's `:soil` layer is total soil water content (mm) in the top
# rooting zone. Multiply by this scale to convert to volumetric water
# content (m^3/m^3) given a 1000 mm rooting depth and 1.3/2.56 bulk-
# density-to-particle-density ratio (NicheMapR convention).
const _TERRACLIMATE_SOIL_TO_VOLUMETRIC = 1.0 / (1000.0 * (1.0 - 1.3 / 2.56))

weather_loader(::Type{<:TerraClimate}) = YearlyTimeSeries()

function weather_variables(::Type{<:TerraClimate})
    (
        WeatherVariable(:maximum_temperature,          :tmax, u"°C"),
        WeatherVariable(:minimum_temperature,          :tmin, u"°C"),
        WeatherVariable(:wind_speed,                   :ws,   u"m/s"),
        WeatherVariable(:vapour_pressure_deficit,      :vpd,  u"kPa"),
        WeatherVariable(:downward_shortwave_radiation, :srad, u"W/m^2"),
        WeatherVariable(:rainfall,                     :ppt,  u"kg/m^2"),
        WeatherVariable(:soil_moisture,                :soil, 1,
                        raw -> raw * _TERRACLIMATE_SOIL_TO_VOLUMETRIC),
    )
end
