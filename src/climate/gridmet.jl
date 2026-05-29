# gridMET (METDATA) bindings — daily ~4 km gridded weather over the
# contiguous United States, 1979-present. One annual NetCDF per layer,
# downloaded and cached via `RasterDataSources.GRIDMET`.
#
# Wind (`:vs`) is reported at 10 m; the chain's
# `derive!(:reference_wind_max)` step applies the same 10 m → 2 m power-law
# shear correction used for NCEP/TerraClimate.
#
# `:rmin` (minimum daily RH, occurring at Tmax) maps to `reference_humidity_min`,
# `:rmax` (maximum daily RH, occurring at Tmin) maps to `reference_humidity_max`.
# A more accurate site-elevation correction would re-evaluate RH against the
# lapse-corrected Tmin/Tmax via vapour-pressure conservation; that's not yet
# wired into the declarative chain.

temporal_resolution(::Type{<:GRIDMET}) = DailyResolution()
weather_loader(::Type{<:GRIDMET}) = YearlyTimeSeries()

function weather_variables(::Type{<:GRIDMET})
    (
        WeatherVariable(:maximum_temperature, :tmmx, u"K"),
        WeatherVariable(:minimum_temperature, :tmmn, u"K"),
        # `:pr` is mm/day depth; numerically equal to kg/m².
        WeatherVariable(:rainfall, :pr, u"kg/m^2"),
        # gridMET RH is stored as percent — convert to 0-1 fraction.
        WeatherVariable(:reference_humidity_min, :rmin, 1, _gridmet_percent_to_fraction),
        WeatherVariable(:reference_humidity_max, :rmax, 1, _gridmet_percent_to_fraction),
        WeatherVariable(:downward_shortwave_radiation, :srad, u"W/m^2"),
        # 10 m wind speed — `derive!(:reference_wind_max)` applies the
        # 10 m → 2 m power-law shear factor.
        WeatherVariable(:wind_speed, :vs, u"m/s"),
    )
end

_gridmet_percent_to_fraction(raw) = raw / 100.0
