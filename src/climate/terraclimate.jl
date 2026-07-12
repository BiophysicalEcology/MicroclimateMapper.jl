# TerraClimate bindings: a yearly NetCDF time series (one file per year, each containing 12 months).

# TerraClimate's `:soil` layer is total soil water content (mm) in the top
# rooting zone. Multiply by this scale to convert to volumetric water
# content (m^3/m^3) given a 1000 mm rooting depth and 1.3/2.56 bulk-
# density-to-particle-density ratio (NicheMapR convention).
_terraclimate_soil_to_volumetric(raw) = raw / (1000.0 * (1.0 - 1.3 / 2.56))

loader(::Type{<:TerraClimate}) = YearlyTimeSeries()

function variables(::Type{<:TerraClimate})
    (
        Variable(Temperature(Maximum()), :tmax, u"°C"),
        Variable(Temperature(Minimum()), :tmin, u"°C"),
        Variable(WindSpeed(), :ws, u"m/s"),
        Variable(VapourPressureDeficit(), :vpd, u"kPa"),
        Variable(GlobalRadiation(), :srad, u"W/m^2"),
        Variable(Rainfall(), :ppt, u"kg/m^2"),
        Variable(SoilMoisture(), :soil, 1, _terraclimate_soil_to_volumetric),
    )
end
