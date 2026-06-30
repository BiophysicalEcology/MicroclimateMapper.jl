# CRUCL2 — CRU CL 2.0 monthly mean climate, 10-minute resolution, 1961–1990.
#
# Single NetCDF `cru_cl2.nc` with a `month` dimension (1–12) for all variables.
# Uses SingleFileBands loader — no per-field file selection needed.
# weather_calendar defaults to Monthly().
#
# The file ships mean temperature, diurnal range and relative humidity rather
# than the extremes and vapour pressure the solver wants. Each is declared as
# the quantity it is; the generic derivation chain then builds temperature_max
# = mean + range/2, temperature_min = mean − range/2, and actual_vapour_pressure
# from relative humidity and mean temperature.

loader(::Type{CRUCL2}) = SingleFileBands()

function variables(::Type{CRUCL2})
    (
        Variable(Temperature(Mean()), :tmp, u"°C"),
        # A diurnal range is a temperature difference: its magnitude is identical
        # in °C and K, so the working unit u"K" carries it with no affine offset.
        Variable(Temperature(DiurnalRange()), :dtr, u"K"),
        # reh is relative humidity percent (0–100); transform → fraction (0–1).
        Variable(RelativeHumidity(), :reh, 1, raw -> raw / 100),
        Variable(WindSpeed(), :wnd, u"m/s"),
        Variable(Rainfall(), :pre, u"kg/m^2"),
        # sunp is sunshine percentage (0–100); transform converts to cloud fraction (0–1).
        Variable(CloudCover(), :sunp, 1, raw -> (100.0 - raw) / 100.0),
        # The source ships its own grid elevation; declared as a quantity it loads
        # through the same path and serves as the lapse-rate reference.
        Variable(Elevation(), :elv, u"m"),
    )
end

# Use the CRUCL2 :elv layer as the DEM. At 10-minute (~18 km) resolution
# terrain is effectively flat, but no lapse rate correction is needed because
# the weather grid and DEM share the same elevation reference.
function _load_dem(::Type{CRUCL2}, area::Extent)
    path     = getraster(CRUCL2)
    buffered = Extents.buffer(area, (X = 0.2, Y = 0.2))   # ≥1 cell at 10-min resolution
    elv      = read(crop(Raster(path; name=:elv, lazy=true); to=buffered, touches=true))
    return Rasters.replace_missing(elv, 0)
end

# load_template for CRUCL2: uses the :elv band from the single NetCDF as the
# run grid. The generic load_template passes extent as an RDS keyword which
# doesn't work for file-path-based sources; dispatch here loads :elv directly.
load_template(T::Type, site::GeocodeResult) = load_template(T, site.extent)
function load_template(::Type{CRUCL2}, extent::Extent)
    path = getraster(CRUCL2)
    return crop(Raster(path; name=:elv, lazy=true); to=extent, touches=true)
end
