# CRUCL2 — CRU CL 2.0 monthly mean climate, 10-minute resolution, 1961–1990.
#
# Single NetCDF `cru_cl2.nc` with a `month` dimension (1–12) for all variables
# plus a static `elv` elevation layer used for lapse rate correction.
# tmax/tmin are derived from tmp ± dtr/2; actual_vapour_pressure from reh × e_sat(tmp);
# cloud_cover from (100 − sunp)/100. temporal_resolution defaults to MonthlyResolution().

# Use the CRUCL2 elv layer as the DEM. At 10-minute (~18 km) resolution
# terrain is effectively flat, but no lapse rate correction is needed because
# the weather grid and DEM share the same elevation reference.
function _load_dem(::Type{CRUCL2}, area::Extent)
    path     = getraster(CRUCL2)
    buffered = Extents.buffer(area, (X = 0.2, Y = 0.2))   # ≥1 cell at 10-min resolution
    elv      = read(crop(Raster(path; name = :elv, lazy = true);
                         to = buffered, touches = true))
    return Rasters.replace_missing(elv, 0)
end

function weather_variables(::Type{CRUCL2})
    (
        WeatherVariable(:maximum_temperature,    :tmax, u"°C"),
        WeatherVariable(:minimum_temperature,    :tmin, u"°C"),
        WeatherVariable(:wind_speed,             :wnd,  u"m/s"),
        WeatherVariable(:rainfall,               :pre,  u"kg/m^2"),
        WeatherVariable(:actual_vapour_pressure, :vapr, u"kPa"),
        WeatherVariable(:cloud_cover,            :cld,  1),
    )
end

# CRUCL2 ships its own grid elevation in the `elv` layer; return it as the
# lapse-rate reference elevation so corrections against the local SRTM DEM work.
weather_grid_elevation(::Type{CRUCL2}, weather, I) =
    Float64(weather[:elv][I..., Ti(1)]) * u"m"

function _load_weather(::Type{CRUCL2}, area::Extent, years)
    path   = getraster(CRUCL2)
    nyears = length(years)

    # Read all monthly variables from the single NetCDF in one pass, then crop.
    raw = read(crop(
        RasterStack(path; name = (:tmp, :dtr, :reh, :sunp, :wnd, :pre), lazy = true);
        to = area, touches = true))

    # Tile (X, Y, 12) → (X, Y, 12*nyears) and rewrap with a Ti axis.
    # CRUCL2 is a fixed 1961–1990 climatology; the same 12 months tile for
    # every requested year. Pattern mirrors MultiBandFutureClimatology in
    # weather.jl and _build_monthly_climatology.
    function tile(layer)
        data = parent(layer)
        size(data, 3) == 12 || error(
            "CRUCL2: expected 12 months in layer, got $(size(data, 3))")
        tiled = nyears == 1 ? data : repeat(data; outer = (1, 1, nyears))
        return Raster(tiled, (dims(layer)[1:2]..., Ti(1:(12 * nyears)));
                      crs = crs(layer))
    end

    tmp  = tile(raw[:tmp])    # mean temperature (°C)
    dtr  = tile(raw[:dtr])    # diurnal temperature range (°C)
    reh  = tile(raw[:reh])    # relative humidity (%)
    sunp = tile(raw[:sunp])   # sunshine percentage (%)
    wnd  = tile(raw[:wnd])    # wind speed (m/s, measured at 10 m)
    pre  = tile(raw[:pre])    # precipitation (mm = kg/m²)
    # TODO: load raw[:wet] (wet days/month) and distribute monthly rainfall
    # across rainy days when splining to a 365-day calendar, matching
    # NicheMapR micro_global.R lines ~799 and ~913–957 (RAINYDAYS / rainfrac).
    # Currently rainfall_daily stays zero in monthly mode.

    sp = dims(tmp)[1:2]
    cr = crs(tmp)
    ti = Ti(1:(12 * nyears))
    wrap(data) = Raster(data, (sp..., ti); crs = cr)

    # tmax and tmin from mean temperature ± half the diurnal range.
    tmax = wrap(parent(tmp) .+ parent(dtr) ./ 2)
    tmin = wrap(parent(tmp) .- parent(dtr) ./ 2)

    # Actual vapour pressure: e_sat(tmp) × reh/100. Use the same GoffGratch
    # equation as the derivation chain so vapr is internally consistent.
    e_sat_kpa = map(parent(tmp)) do t_c
        ustrip(u"kPa", vapour_pressure(GoffGratch(), (t_c + 273.15) * u"K"))
    end
    vapr = wrap(e_sat_kpa .* parent(reh) ./ 100.0)

    # Cloud cover from sunshine percentage: cld = (100 - sunp) / 100.
    cld = wrap((100.0 .- parent(sunp)) ./ 100.0)

    # Load elevation separately (no month dimension) and tile to Ti so the
    # returned stack has uniform dimensions. The `weather_grid_elevation`
    # dispatch reads :elv[I..., Ti(1)] for lapse rate correction.
    elv_raw  = read(crop(Raster(path; name = :elv, lazy = true);
                         to = area, touches = true))
    elv_data = repeat(parent(elv_raw); outer = (1, 1, 12 * nyears))
    elv      = Raster(elv_data, (dims(elv_raw)..., ti); crs = crs(elv_raw))

    stack = RasterStack((; tmax, tmin, wnd, pre, vapr, cld, elv))
    return Rasters.replace_missing(stack, NaN)
end
