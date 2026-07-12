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
# the actual-vapour-pressure derivation is skipped automatically because
# relative humidity is declared as a native. Pressure comes from mean sea
# level pressure (`:psl`) — declared here as `Pressure()`, so no explicit
# sea-level-to-site pressure correction is applied yet. That derivation
# should be added when the full BARRA end-to-end path lands.

weather_calendar(::Type{<:BARRA}) = Daily()
native_timestep(::Type{<:BARRA}) = Hourly()
loader(::Type{<:BARRA}) = MonthlyTimeSeries()

function variables(::Type{<:BARRA})
    (
        Variable(Reference(Temperature()), :tas, u"K"),
        Variable(WindSpeed(), :sfcWind, u"m/s"),
        Variable(Reference(RelativeHumidity()), :hurs, 1, _barra_percent_to_fraction),
        # `:psl` is mean sea level pressure — see file header note; treated as
        # site pressure here pending an explicit pressure-from-sea-level
        # correction derivation.
        Variable(Pressure(), :psl, u"Pa"),
        Variable(GlobalRadiation(), :rsds, u"W/m^2"),
        Variable(LongwaveRadiation(), :rlds, u"W/m^2"),
        # `:pr` (kg/m²/s) is the Day-frequency mean rate; × 86400 s → the
        # day's total. Currently loaded at hourly frequency alongside the
        # rest; the "single midnight event" placement from main is TODO.
        Variable(Rainfall(), :pr, u"kg/m^2", raw -> raw * 86400.0),
        # BARRA's own :orog grid — declared as a static Elevation, so
        # weather_grid_elevation picks it up automatically.
        Variable(Elevation(), :orog, u"m"),
    )
end

_barra_percent_to_fraction(raw) = raw / 100.0

# BARRA's own `:orog` grid as the DEM — same role as CRUCL2's `:elv`. 0.3°
# buffer matches SRTM's, comfortably wider than a BARRA-R2 cell (~11 km).
function _load_dem(T::Type{<:BARRA}, area::Extent)
    buffered = Extents.buffer(area, (X = 0.3, Y = 0.3))
    path = getraster(T, :orog)
    orog = read(crop(Raster(path; name = :orog, lazy = true); to = buffered, touches = true))
    return Rasters.replace_missing(orog, 0)
end
