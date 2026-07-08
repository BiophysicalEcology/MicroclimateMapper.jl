# CRUCL2 — CRU CL 2.0 monthly mean climate, 10-minute resolution, 1961–1990.
#
# Single NetCDF `cru_cl2.nc` with a `month` dimension (1–12) for all variables.
# Uses SingleFileBands loader — no per-field file selection needed.
# tmax/tmin and vapr require two-input derivations (WeatherVariable.transform
# only handles one field), so they are computed in _post_load_stack!.
# temporal_resolution defaults to MonthlyResolution().

weather_loader(::Type{CRUCL2}) = SingleFileBands()

# Primary fields to load from the single NetCDF. :tmp, :dtr, :reh are
# intermediate inputs for the two-input derivations in _post_load_stack!
# and do not appear in weather_variables. :elv tiles with the rest so the
# returned stack has uniform dimensions for weather_grid_elevation.
primary_layers(::Type{CRUCL2}) = (:tmp, :dtr, :reh, :wnd, :pre, :sunp, :elv)

function weather_variables(::Type{CRUCL2})
    (
        WeatherVariable(:maximum_temperature,    :tmax, u"°C"),
        WeatherVariable(:minimum_temperature,    :tmin, u"°C"),
        WeatherVariable(:wind_speed,             :wnd,  u"m/s"),
        WeatherVariable(:rainfall,               :pre,  u"kg/m^2"),
        WeatherVariable(:actual_vapour_pressure, :vapr, u"kPa"),
        # sunp is sunshine percentage (0–100); transform converts to cloud fraction (0–1).
        WeatherVariable(:cloud_cover, :sunp, 1, s -> (100.0 - s) / 100.0),
    )
end

# CRUCL2 ships its own grid elevation in the :elv layer; return it as the
# lapse-rate reference elevation so corrections against the local SRTM DEM work.
weather_grid_elevation(::Type{CRUCL2}, weather, I) =
    Float64(weather[:elv][I..., Ti(1)]) * u"m"

# Compute tmax, tmin, vapr from the raw intermediate fields loaded by
# SingleFileBands, then rebuild the stack with only the canonical
# native fields (removing :tmp, :dtr, :reh).
function _post_load_stack!(::Type{CRUCL2}, stack, _)
    ti = dims(stack[:tmp], Ti)
    other = otherdims(stack[:tmp], Ti)
    cr = crs(stack[:tmp])
    wrap(data) = Raster(data, (other..., ti); crs = cr)

    tmp = parent(stack[:tmp])
    dtr = parent(stack[:dtr])
    reh = parent(stack[:reh])

    tmax = wrap(tmp .+ dtr ./ 2)
    tmin = wrap(tmp .- dtr ./ 2)
    # Actual vapour pressure: e_sat(tmp) × reh/100.
    # GoffGratch matches the derivation chain so vapr is internally consistent.
    vapr = wrap(map(t -> ustrip(u"kPa", vapour_pressure(GoffGratch(), (t + 273.15) * u"K")), tmp)
                .* reh ./ 100)

    return RasterStack((; tmax, tmin,
                          wnd  = stack[:wnd],
                          pre  = stack[:pre],
                          vapr,
                          sunp = stack[:sunp],
                          elv  = stack[:elv]))
end

# Use the CRUCL2 :elv layer as the DEM. At 10-minute (~18 km) resolution
# terrain is effectively flat, but no lapse rate correction is needed because
# the weather grid and DEM share the same elevation reference.
function _load_dem(::Type{CRUCL2}, area::Extent)
    path     = getraster(CRUCL2)
    buffered = Extents.buffer(area, (X = 0.2, Y = 0.2))   # ≥1 cell at 10-min resolution
    elv      = read(crop(Raster(path; name = :elv, lazy = true);
                         to = buffered, touches = true))
    return Rasters.replace_missing(elv, 0)
end

# load_template for CRUCL2: uses the :elv band from the single NetCDF as the
# run grid. The generic load_template passes extent as an RDS keyword which
# doesn't work for file-path-based sources; dispatch here loads :elv directly.
load_template(::Type{CRUCL2}, site::GeocodeResult) = load_template(CRUCL2, site.extent)
function load_template(::Type{CRUCL2}, extent::Extent)
    path = getraster(CRUCL2)
    read(crop(Raster(path; name = :elv, lazy = true); to = extent, touches = true))
end
