# gridMET (METDATA) bindings — daily ~4 km gridded weather over the

weather_calendar(::Type{<:GRIDMET}) = Daily()
loader(::Type{<:GRIDMET}) = YearlyTimeSeries()

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
        # 10 m wind speed — `derive!(:reference_wind_max)` applies the
        Variable(WindSpeed(), :vs, u"m/s"),
    )
end

# NetCDF variable names in gridMET's files
file_variable(::Type{<:GRIDMET}, name::Symbol) = get((
    tmmx = :air_temperature, tmmn = :air_temperature, pr = :precipitation_amount,
    rmin = :relative_humidity, rmax = :relative_humidity, srad = :surface_downwelling_shortwave_flux_in_air,
    vs = :wind_speed, elev = :elevation), name, name)
