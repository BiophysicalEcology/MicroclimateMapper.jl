# BARRA bindings — Bureau of Meteorology's high-resolution Australian
# regional reanalysis (BARRA-R2 ~11 km / BARRA-C2 ~4 km), Australia's
# analog to ERA5. One NetCDF per variable per real calendar month.
#
# Hourly met variables are loaded at native `Hour` frequency (like ERA5).
# Rainfall instead comes from the `Day`-frequency file and is placed as a
# single event at each day's first hour (midnight), matching the "one event
# at midnight" convention used elsewhere for daily-total rainfall. BARRA's
# own `:orog` static elevation grid backs `weather_grid_elevation`, same
# role as CRUCL2's bundled `:elv` and GRIDMET's `:elev`.
#
# BARRA gives relative humidity (`:hurs`) directly rather than dewpoint, so
# it uses `_DERIVATIONS_HOURLY_NATIVE_RH` instead of ERA5's dewpoint-based
# chain. Pressure comes from mean sea level pressure (`:psl`) rather than
# native surface pressure, adjusted to the DEM site elevation via
# `derive!(:pressure_from_sea_level)` — BARRA's own grid elevation
# generally differs from the simulated site's.

temporal_resolution(::Type{<:BARRA}) = HourlyResolution()
weather_derivations(::Type{<:BARRA}) = (_DERIVATIONS_HOURLY_NATIVE_RH..., Val(:pressure_from_sea_level))

function weather_variables(::Type{<:BARRA})
    (
        WeatherVariable(:reference_temperature, :tas, u"K"),
        WeatherVariable(:wind_speed, :sfcWind, u"m/s"),
        WeatherVariable(:reference_humidity, :hurs, 1, _barra_percent_to_fraction),
        WeatherVariable(:sea_level_pressure, :psl, u"Pa"),
        WeatherVariable(:global_radiation, :rsds, u"W/m^2"),
        WeatherVariable(:longwave_radiation, :rlds, u"W/m^2"),
        # `:pr` (kg/m²/s) is the Day-frequency mean rate; × 86400 s → the
        # day's total, placed at midnight by `_load_weather` below.
        WeatherVariable(:rainfall, :pr, u"kg/m^2", raw -> raw * 86400.0),
    )
end

_barra_percent_to_fraction(raw) = raw / 100.0

weather_grid_elevation(::Type{<:BARRA}, weather, I) =
    Float64(weather[:orog][I..., Ti(1)]) * u"m"

# BARRA's own `:orog` grid as the DEM — same role as CRUCL2's `:elv`. 0.3°
# buffer matches SRTM's, comfortably wider than a BARRA-R2 cell (~11 km).
function _load_dem(T::Type{<:BARRA}, area::Extent)
    buffered = Extents.buffer(area, (X = 0.3, Y = 0.3))
    path = getraster(T, :orog)
    orog = read(crop(Raster(path; name = :orog, lazy = true); to = buffered, touches = true))
    return Rasters.replace_missing(orog, 0)
end

function _load_weather(T::Type{<:BARRA{P, D}}, area::Extent, years) where {P, D}
    met = _load_layers(MonthlyTimeSeries(), T,
        (:tas, :sfcWind, :hurs, :psl, :rsds, :rlds), area, years)
    ref = first(values(met))

    daily = _load_layers(MonthlyTimeSeries(), BARRA{P, D, Day}, (:pr,), area, years)
    pr = _daily_to_hourly_midnight(daily.pr, ref)

    orog_path = getraster(T, :orog)
    orog_raw = read(crop(Raster(orog_path; name = :orog, lazy = true); to = area, touches = true))
    orog = _static_to_ti(orog_raw, ref)

    stack = RasterStack(merge(met, (; pr, orog)))
    return Rasters.replace_missing(stack, NaN)
end

supports_points_loading(::Type{<:BARRA}) = true

# Points-native counterpart of `_load_weather`, used by `MicroVectorProblem`.
function _load_weather_points(T::Type{<:BARRA{P, D}}, points_dim, years) where {P, D}
    met = _load_layers_at_points(MonthlyTimeSeries(), T,
        (:tas, :sfcWind, :hurs, :psl, :rsds, :rlds), points_dim, years)
    ref = first(values(met))

    daily = _load_layers_at_points(MonthlyTimeSeries(), BARRA{P, D, Day}, (:pr,), points_dim, years)
    pr = _daily_to_hourly_midnight_points(daily.pr, ref)

    orog = _static_to_ti_points(_load_field_at_points(T, :orog, points_dim), ref)

    stack = RasterStack(merge(met, (; pr, orog)))
    return Rasters.replace_missing(stack, NaN)
end
