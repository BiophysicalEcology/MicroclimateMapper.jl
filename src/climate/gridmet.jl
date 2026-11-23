# gridMET (METDATA) bindings — daily ~4 km gridded weather over the
# contiguous US. File variable names come from `RasterDataSources.layerkeys`.

weather_calendar(::Type{<:GRIDMET}) = Daily()
loader(::Type{<:GRIDMET}) = YearlyTimeSeries()
# gridMET's elevation grid is its own source.
static_source(::Type{GRIDMET}) = GRIDMET{RasterDataSources.Elevation}

function variables(::Type{<:GRIDMET})
    (
        Variable(Temperature(Maximum()), :tmmx, u"K"),
        Variable(Temperature(Minimum()), :tmmn, u"K"),
        # `:pr` is mm/day depth; numerically equal to kg/m².
        Variable(Rainfall(), :pr, u"kg/m^2"),
        # gridMET relative humidity is stored as percent — convert to 0-1 fraction.
        Variable(Reference(RelativeHumidity(Minimum())), :rmin, 1, percent_to_fraction),
        Variable(Reference(RelativeHumidity(Maximum())), :rmax, 1, percent_to_fraction),
        Variable(GlobalRadiation(), :srad, u"W/m^2"),
        # 10 m wind speed — `derive!(:reference_wind_max)` applies the height correction.
        Variable(WindSpeed(), :vs, u"m/s"),
        Variable(Elevation(), :elev, u"m"),
    )
end
