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

_to_float_with_nan(raster::Raster) =
    map(x -> ismissing(x) ? NaN : Float64(x), parent(raster))

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
function compute_terrain_grids(utm_dem::Raster; n_horizon_angles = 24)
    canonical_dem  = _canonicalize(utm_dem)
    canonical_dims = dims(canonical_dem)
    raster_crs     = crs(canonical_dem)

    elevation_matrix = _to_float_with_nan(canonical_dem)
    x_coordinates    = collect(lookup(canonical_dem, X))
    y_coordinates    = collect(lookup(canonical_dem, Y))
    cell_size        = (abs(x_coordinates[2] - x_coordinates[1]),
                        abs(y_coordinates[2] - y_coordinates[1]))

    slope_matrix  = Geomorphometry.slope( elevation_matrix; method = Horn(), cellsize = cell_size)
    aspect_matrix = Geomorphometry.aspect(elevation_matrix; method = Horn(), cellsize = cell_size)

    horizon_matrix = _compute_horizon_angles_matrix(
        elevation_matrix, x_coordinates, y_coordinates, n_horizon_angles,
    )

    latlon_dem  = Rasters.resample(canonical_dem; crs = GeoFormatTypes.EPSG(4326))
    lat_south, lat_north = extrema(lookup(latlon_dem, Y))
    lon_west,  lon_east  = extrema(lookup(latlon_dem, X))
    latitudes  = collect(range(lat_south, lat_north; length = length(y_coordinates)))
    longitudes = collect(range(lon_west,  lon_east;  length = length(x_coordinates)))
    latitude_matrix  = repeat(reshape(latitudes,  1, :), length(x_coordinates), 1)  # (nx, ny)
    longitude_matrix = repeat(longitudes, 1, length(y_coordinates))                  # (nx, ny)

    tag_2d(matrix, unit) =
        Raster(map(x -> isnan(x) ? missing : x * unit, matrix), canonical_dims; crs = raster_crs)

    elevation_m   = tag_2d(elevation_matrix, 1.0u"m")
    slope_deg     = tag_2d(slope_matrix,     1.0u"°")
    aspect_deg    = tag_2d(aspect_matrix,    1.0u"°")
    latitude_deg  = tag_2d(latitude_matrix,  1.0u"°")
    longitude_deg = tag_2d(longitude_matrix, 1.0u"°")
    pressure_pa   = map(e -> ismissing(e) ? missing : atmospheric_pressure(e), elevation_m)

    horizon_dims = (canonical_dims..., Dim{:azimuth}(1:n_horizon_angles))
    horizon_angles_deg = Raster(
        map(x -> isnan(x) ? missing : x * u"°", horizon_matrix),
        horizon_dims;
        crs = raster_crs,
    )

    return (; elevation_m, slope_deg, aspect_deg,
              latitude_deg, longitude_deg, pressure_pa,
              horizon_angles_deg)
end

"""
    compute_horizon_angles(elevation::Raster, n_directions; verbose=true) -> Raster

Compute terrain horizon-elevation angles for `n_directions` evenly-spaced
azimuth directions. Returns a 3-D `Raster` with dimensions
`(X, Y, Dim{:azimuth})` and unit `u"°"`.

`elevation` must be a `Raster` carrying an elevation field (with or without
units); missing/`NaN` cells are propagated as `missing`.
"""
function compute_horizon_angles(elevation::Raster, n_directions; verbose = true)
    canonical_elevation = _canonicalize(elevation)
    elevation_matrix    = _to_float_with_nan(canonical_elevation)
    x_coordinates       = collect(lookup(canonical_elevation, X))
    y_coordinates       = collect(lookup(canonical_elevation, Y))

    horizon_matrix = _compute_horizon_angles_matrix(
        elevation_matrix, x_coordinates, y_coordinates, n_directions; verbose,
    )

    return Raster(
        map(x -> isnan(x) ? missing : x * u"°", horizon_matrix),
        (dims(canonical_elevation)..., Dim{:azimuth}(1:n_directions));
        crs = crs(canonical_elevation),
    )
end

# Kernel: operates on a plain `(nx, ny)` Float64 matrix with both coordinate
# axes ascending. Returns a `(nx, ny, n_directions)` Float64 array of horizon
# elevation angles in degrees (NaN where elevation is NaN).
function _compute_horizon_angles_matrix(
    elevation_matrix::AbstractMatrix{Float64},
    x_coordinates::AbstractVector,
    y_coordinates::AbstractVector,
    n_directions::Integer;
    verbose = true,
)
    nx, ny  = size(elevation_matrix)
    horizons = zeros(Float64, nx, ny, n_directions)
    dx = Float64(x_coordinates[2] - x_coordinates[1])
    dy = Float64(y_coordinates[2] - y_coordinates[1])

    for direction in 1:n_directions
        azimuth      = 2π * (direction - 1) / n_directions
        raw_step_x   = sin(azimuth) / abs(dx)
        raw_step_y   = cos(azimuth) / dy
        normaliser   = max(abs(raw_step_x), abs(raw_step_y))
        normaliser ≈ 0 && continue
        step_x = raw_step_x / normaliser
        step_y = raw_step_y / normaliser

        verbose && print("  horizon direction $direction/$n_directions  \r")
        for j in 1:ny, i in 1:nx
            if isnan(elevation_matrix[i, j])
                horizons[i, j, direction] = NaN
                continue
            end
            elevation_here = elevation_matrix[i, j]
            max_tan = 0.0
            step    = 1
            while true
                xi = round(Int, i + step * step_x)
                yi = round(Int, j + step * step_y)
                (xi < 1 || xi > nx || yi < 1 || yi > ny) && break
                if !isnan(elevation_matrix[xi, yi])
                    distance = sqrt((x_coordinates[xi] - x_coordinates[i])^2 +
                                    (y_coordinates[yi] - y_coordinates[j])^2)
                    if distance > 0
                        tangent = (elevation_matrix[xi, yi] - elevation_here) / distance
                        tangent > max_tan && (max_tan = tangent)
                    end
                end
                step += 1
            end
            horizons[i, j, direction] = atand(max_tan)
        end
    end
    verbose && println()
    return horizons
end
