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

# :elev is a static file on its own grid — resample onto the loaded layers'
# grid and tile across Ti (as CRUCL2 does for :elv) so it slices cleanly
# alongside the other layers.
function _post_load_stack!(::Type{<:GRIDMET}, stack, _years)
    template = first(values(stack))
    path = getraster(GRIDMET, :elev)
    elev_raw = read(Raster(path; lazy = true)[Dim{:day}(1)])
    elev_2d = Rasters.replace_missing(
        Rasters.resample(elev_raw; to = dims(template, (X, Y))), NaN)
    ti = dims(template, Ti)
    tiled = repeat(reshape(parent(elev_2d), size(elev_2d)..., 1); outer = (1, 1, length(ti)))
    elev = Raster(tiled, (dims(elev_2d)..., ti); crs = crs(elev_2d))
    names = keys(stack)
    return RasterStack(NamedTuple{(names..., :elev)}((map(n -> stack[n], names)..., elev)))
end

weather_grid_elevation(::Type{<:GRIDMET}, weather, I) =
    Float64(weather[:elev][I..., Ti(1)]) * u"m"

# Points-mode counterpart of `_post_load_stack!` -- extracts :elev directly at
# each point instead of resampling a separate grid (no X/Y dims to resample onto).
function _post_load_stack_points!(::Type{<:GRIDMET}, stack, _years, points_dim)
    template = first(values(stack))
    ti = dims(template, Ti)
    path = getraster(GRIDMET, :elev)
    elev_pts = _extract_lazy_at_points(Raster(path; lazy = true)[Dim{:day}(1)], points_dim)
    tiled = repeat(reshape(parent(elev_pts), :, 1); outer = (1, length(ti)))
    elev = Raster(tiled, (points_dim, ti))
    names = keys(stack)
    return RasterStack(NamedTuple{(names..., :elev)}((map(n -> stack[n], names)..., elev)))
end
