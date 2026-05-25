"""
    _load_dem(::Type{<:RasterDataSources.RasterDataSource}, area::Extent) -> Raster

Load a DEM for the given lon/lat `area`. Source-specific methods live next
to this declaration (e.g. `terrain/srtm.jl` for `SRTM`).
"""
function _load_dem end

"""
    compute_terrain_grids(utm_dem::Raster; n_horizon_angles=24) -> NamedTuple

Compute the full set of per-pixel terrain rasters needed for solar radiation
and microclimate simulations from a UTM-projected DEM raster.

All returned grids are `Raster`s sharing the dimensions of `utm_dem`
(canonicalised to `(X, Y)` order with ascending lookups). Pixel access is by
named dimension; callers never need to know the storage layout.

# Returns
RasterStack with fields (all `Raster`s):
- `elevation_m`     : 2-D elevation              (`u"m"`)
- `slope_deg`       : 2-D slope                  (`u"°"`)
- `aspect_deg`      : 2-D aspect                 (`u"°"`)
- `latitude_deg`    : 2-D per-pixel latitude     (`u"°"`)
- `longitude_deg`   : 2-D per-pixel longitude    (`u"°"`)
- `pressure_pa`     : 2-D atmospheric pressure   (`u"Pa"`)
- `horizon_angles_deg` : 3-D horizon-elevation angles indexed by
                         `(X, Y, Dim{:azimuth})` (`u"°"`)
"""
function compute_terrain_grids(dem::Raster; n_horizon_angles = 24)
    # Treat the encoded no-data as sea level (0) and drop the missingval —
    # downstream operations (`.* u"m"`, slope/aspect, RasterStack) then have
    # nothing Missing-shaped to reason about.
    # TODO: this is dodgy
    canonical_dem = rebuild(_canonicalize(dem); missingval=nothing)
    # Geomorphometry's Rasters extension returns `(dx, dy)` in metres for
    # both projected and geographic CRSes (converting deg→m at the centre
    # latitude in the lat/lon case). Note: NOT `Rasters.cellarea`, which
    # returns *area* (m²) — wrong shape for the slope/aspect kernels.
    cellsize = Geomorphometry.cellsize(canonical_dem)

    # Slope/aspect stay as Rasters — Geomorphometry's `similar(dem, Float32)`
    # carries the Raster wrapper through.
    slope_raster  = Geomorphometry.slope( canonical_dem; method=Horn(), cellsize)
    aspect_raster = Geomorphometry.aspect(canonical_dem; method=Horn(), cellsize)

    horizon_angles = compute_horizon_angles(canonical_dem; directions = n_horizon_angles)

    # Lat/lon are just the WGS84 lookups broadcast over the grid, tagged
    # with `u"°"` so they match `SolarRadiation.solar_geometry`'s
    # `latitude::Quantity` signature downstream. `DimPoints` yields
    # `(X, Y)` per pixel — X is longitude, Y is latitude.
    raster_dims = dims(canonical_dem)
    raster_crs  = crs(canonical_dem)
    longitude = first.(DimPoints(canonical_dem)) .* u"°"
    latitude  = last.(DimPoints(canonical_dem))  .* u"°"

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
        directions, cellsize = Geomorphometry.cellsize(canonical))
    # Pack the trailing direction axis into a per-pixel SVector via reinterpret.
    permuted = permutedims(raw, (3, 1, 2)) .* u"°"
    data = reinterpret(reshape, SVector{directions, eltype(permuted)}, permuted)
    return Raster(collect(data), dims(canonical); crs = crs(canonical))
end

# Isotropic sky-view factor (Dozier & Frew 1990, flat-surface limit):
# V_sky = (1/N) Σ cos²(h_i). `V_sky = 1` for an unobstructed horizon.
@inline _sky_view_from_horizon(horizon_angles) =
    sum(h -> cos(h)^2, horizon_angles) / length(horizon_angles)

# Terrain utility functions for biophysical grid analyses.
#
# All public functions take and return `Rasters.Raster`s, so callers index
# pixels by named dimension (`r[X=i, Y=j]`) and never by raw storage order.

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
