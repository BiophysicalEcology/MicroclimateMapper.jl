# Terrain utility functions for biophysical grid analyses.
#
# All public functions take and return `Rasters.Raster`s, so callers index
# pixels by named dimension (`r[X=i, Y=j]`) and never by raw storage order.

"""
    get_utm_crs(raster) -> GeoFormatTypes.EPSG

Return the UTM EPSG code whose zone contains the centre of `raster`.
Northern-hemisphere zones are EPSG 326xx; southern-hemisphere zones 327xx.
"""
function get_utm_crs(raster)
    center_lon = lookup(raster, X)[size(raster, X) ÷ 2]
    center_lat = lookup(raster, Y)[size(raster, Y) ÷ 2]
    utm_zone   = floor(Int, (center_lon + 180) / 6) + 1
    epsg_code  = center_lat ≥ 0 ? "326$(lpad(utm_zone, 2, '0'))" :
                                   "327$(lpad(utm_zone, 2, '0'))"
    return GeoFormatTypes.EPSG(parse(Int, epsg_code))
end

# Canonicalise a Raster to `(X, Y)` dim order with both axes ascending.
# Internal helper used by the terrain algorithms that need a fixed orientation.
function _canonicalize(r::Raster)
    r2 = dimnum(r, X) > dimnum(r, Y) ? permutedims(r, (X, Y)) : r
    y_lookup = lookup(r2, Y)
    if length(y_lookup) > 1 && first(y_lookup) > last(y_lookup)
        r2 = reverse(r2; dims = Y)
    end
    return r2
end


"""
    load_utm_dem(center_lon, center_lat, extent_lon, extent_lat) -> NamedTuple

Download an SRTM DEM tile, crop to the requested bounding box, and reproject
to the local UTM zone.

# Returns
NamedTuple with fields:
- `utm_dem`       : reprojected `Raster`
- `x_coords_utm`  : easting vector (m)
- `y_coords_utm`  : northing vector (m)
- `nx_utm`        : number of columns
- `ny_utm`        : number of rows
- `cs`            : `(dx, dy)` cell size tuple (m)
"""
function load_utm_dem(center_lon, center_lat, extent_lon, extent_lat)
    lon_min = center_lon - extent_lon / 2
    lon_max = center_lon + extent_lon / 2
    lat_min = center_lat - extent_lat / 2
    lat_max = center_lat + extent_lat / 2

    tile_paths  = getraster(SRTM; bounds = (lon_min, lat_min, lon_max, lat_max))
    valid_paths = filter(!ismissing, vec(tile_paths))
    isempty(valid_paths) && error("No SRTM tile found for bounds " *
        "($(lon_min), $(lat_min), $(lon_max), $(lat_max)).")
    dem_full  = Raster(only(valid_paths))
    dem_wgs84 = dem_full[X(Between(lon_min, lon_max)), Y(Between(lat_min, lat_max))]

    utm_crs      = get_utm_crs(dem_wgs84)
    utm_dem      = Rasters.resample(dem_wgs84; crs = utm_crs, method = :bilinear)
    x_coords_utm = collect(lookup(utm_dem, X))
    y_coords_utm = collect(lookup(utm_dem, Y))
    nx_utm       = length(x_coords_utm)
    ny_utm       = length(y_coords_utm)
    cs           = (abs(x_coords_utm[2] - x_coords_utm[1]),
                    abs(y_coords_utm[2] - y_coords_utm[1]))

    return (; utm_dem, x_coords_utm, y_coords_utm, nx_utm, ny_utm, cs)
end

"""
    compute_terrain_grids(utm_dem::Raster; n_horizon_angles=24) -> NamedTuple

Compute the full set of per-pixel terrain rasters needed for solar radiation
and microclimate simulations from a UTM-projected DEM raster.

All returned grids are `Raster`s sharing the dimensions of `utm_dem`
(canonicalised to `(X, Y)` order with ascending lookups). Pixel access is by
named dimension; callers never need to know the storage layout.

# Returns
NamedTuple with fields (all `Raster`s):
- `elevation_m`     : 2-D elevation              (`u"m"`)
- `slope_deg`       : 2-D slope                  (`u"°"`)
- `aspect_deg`      : 2-D aspect                 (`u"°"`)
- `latitude_deg`    : 2-D per-pixel latitude     (`u"°"`)
- `longitude_deg`   : 2-D per-pixel longitude    (`u"°"`)
- `pressure_pa`     : 2-D atmospheric pressure   (`u"Pa"`)
- `horizon_angles_deg` : 3-D horizon-elevation angles indexed by
                         `(X, Y, Dim{:azimuth})` (`u"°"`)
"""
const _EARTH_RADIUS_M = 6_378_137.0  # WGS84 equatorial radius

"""
Per-pixel metric cell size for a WGS84 raster, using the mid-latitude as the
representative latitude for the whole grid. Approximate; fine for areas a few
degrees across and away from the poles.
"""
function _wgs84_cellsize(dem::Raster)
    xs = lookup(dem, X)
    ys = lookup(dem, Y)
    dlon_deg = abs(xs[2] - xs[1])
    dlat_deg = abs(ys[2] - ys[1])
    midlat   = (first(ys) + last(ys)) / 2
    dy_m     = _EARTH_RADIUS_M * deg2rad(dlat_deg)
    dx_m     = _EARTH_RADIUS_M * deg2rad(dlon_deg) * cosd(midlat)
    return (dx_m, dy_m)
end

function compute_terrain_grids(dem::Raster; n_horizon_angles = 24)
    # Treat the encoded no-data as sea level (0) and drop the missingval —
    # downstream operations (`.* u"m"`, slope/aspect, RasterStack) then have
    # nothing Missing-shaped to reason about.
    canonical_dem = rebuild(_canonicalize(dem); missingval = nothing)
    cell_size     = _wgs84_cellsize(canonical_dem)

    # Slope/aspect stay as Rasters — Geomorphometry's `similar(dem, Float32)`
    # carries the Raster wrapper through.
    slope_raster  = Geomorphometry.slope( canonical_dem; method = Horn(), cellsize = cell_size)
    aspect_raster = Geomorphometry.aspect(canonical_dem; method = Horn(), cellsize = cell_size)

    horizon_angles = compute_horizon_angles(canonical_dem; directions = n_horizon_angles)

    # Lat/lon are just the WGS84 lookups broadcast over the grid.
    xs = lookup(canonical_dem, X)
    ys = lookup(canonical_dem, Y)
    nx, ny = length(xs), length(ys)
    raster_dims = dims(canonical_dem)
    raster_crs  = crs(canonical_dem)
    latitude  = Raster(repeat(reshape(collect(ys), 1, ny), nx, 1) .* u"°",
                       raster_dims; crs = raster_crs)
    longitude = Raster(repeat(collect(xs), 1, ny) .* u"°",
                       raster_dims; crs = raster_crs)

    # Layer names match `Site` field names so a per-pixel stack view splats
    # directly into the Site constructor.
    elevation = canonical_dem .* u"m"
    slope     = slope_raster  .* u"°"
    aspect    = aspect_raster .* u"°"
    pressure  = atmospheric_pressure.(elevation)

    return RasterStack((;
        elevation, slope, aspect, latitude, longitude,
        atmospheric_pressure = pressure, horizon_angles,
    ))
end

"""
    compute_horizon_angles(elevation::Raster; directions=16) -> Raster

Wraps `Geomorphometry.horizon_angle` to return a 2-D `Raster` whose values
are length-`directions` `SVector`s of horizon angles (in `u"°"`).
"""
function compute_horizon_angles(elevation::Raster; directions = 16)
    canonical = _canonicalize(elevation)
    raw = Geomorphometry.horizon_angle(parent(canonical);
        directions, cellsize = _wgs84_cellsize(canonical))
    # Pack the trailing direction axis into a per-pixel SVector via reinterpret.
    permuted = permutedims(raw, (3, 1, 2)) .* u"°"
    data = reinterpret(reshape, SVector{directions, eltype(permuted)}, permuted)
    return Raster(collect(data), dims(canonical); crs = crs(canonical))
end

# Isotropic sky-view factor (Dozier & Frew 1990, flat-surface limit):
# V_sky = (1/N) Σ cos²(h_i). `V_sky = 1` for an unobstructed horizon.
@inline _sky_view_from_horizon(horizon_angles) =
    sum(h -> cos(h)^2, horizon_angles) / length(horizon_angles)
