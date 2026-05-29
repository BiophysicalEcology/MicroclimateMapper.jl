"""
    _load_dem(::Type{<:RasterDataSources.RasterDataSource}, area::Extent) -> Raster

Load a DEM for the given lon/lat `area`. Source-specific methods live next
to this declaration (e.g. `terrain/srtm.jl` for `SRTM`).
"""
function _load_dem end

"""
    compute_terrain_grids(dem::Raster;
                          template = nothing,
                          n_horizon_angles = 16) -> RasterStack

Compute the full set of per-pixel terrain rasters needed for solar radiation
and microclimate simulations.

Slope, aspect, and horizon-angle vectors are always computed on `dem` at
its **native** resolution — coarsening first would smear out the relief
that drives shadow / slope-radiation effects. When `template` is given
each result is resampled onto the template grid (`:average` for
elevation / slope / per-azimuth horizon; sin/cos circular mean for
aspect), so the returned stack lives at the coarse run resolution while
still reflecting fine-scale terrain.

# Returns
`RasterStack` with layers (all `Raster`s):
- `elevation`            : 2-D elevation              (`u"m"`)
- `slope`                : 2-D slope                  (`u"°"`)
- `aspect`               : 2-D aspect                 (`u"°"`)
- `latitude`             : 2-D per-pixel latitude     (`u"°"`)
- `longitude`            : 2-D per-pixel longitude    (`u"°"`)
- `atmospheric_pressure` : 2-D atmospheric pressure   (`u"Pa"`)
- `horizon_angles`       : 2-D `SVector{n_horizon_angles, °}` per pixel
"""
function compute_terrain_grids(dem::Raster;
                               template = nothing,
                               n_horizon_angles = 16)
    canonical_native = rebuild(_canonicalize(dem); missingval = nothing)
    # `Geomorphometry.cellsize` rejects Points sampling — DEMs sometimes
    # arrive with Points and need to be flipped to Intervals(Center()).
    canonical_native = Rasters.set(canonical_native,
                                   X => Intervals(Center()),
                                   Y => Intervals(Center()))
    cellsize_native = Geomorphometry.cellsize(canonical_native)

    slope_native = Geomorphometry.slope( canonical_native; method = Horn(), cellsize = cellsize_native)
    aspect_native = Geomorphometry.aspect(canonical_native; method = Horn(), cellsize = cellsize_native)
    horizon_native_stack = _compute_horizon_native(canonical_native; directions = n_horizon_angles)

    target = template === nothing ? canonical_native : template

    elevation = Rasters.resample(canonical_native; to = target, method = :average) .* u"m"
    slope = Rasters.resample(slope_native; to = target, method = :average) .* u"°"
    aspect = _circular_mean_aspect(aspect_native, target)                       .* u"°"
    horizon = _resample_horizon(horizon_native_stack, target, n_horizon_angles)

    longitude = first.(DimPoints(elevation)) .* u"°"
    latitude = last.(DimPoints(elevation))  .* u"°"
    pressure = atmospheric_pressure.(elevation)

    return RasterStack((;
        elevation, slope, aspect, latitude, longitude,
        atmospheric_pressure = pressure, horizon_angles = horizon,
    ))
end

# Aspect is circular (0° ≡ 360°) in compass convention (clockwise from
# north, [0, 360)). Plain mean would average 350°/10° as 180° (south);
# instead average sin/cos and atan2 back, then re-wrap into [0, 360) so
# downstream solar geometry sees the same convention Geomorphometry
# emits at native resolution.
function _circular_mean_aspect(aspect_native::Raster, target)
    target === aspect_native && return aspect_native
    asin = Rasters.resample(sind.(aspect_native); to = target, method = :average)
    acos = Rasters.resample(cosd.(aspect_native); to = target, method = :average)
    return mod.(atand.(asin, acos), 360)
end

# Compute horizon angles at native resolution and return them as a (X, Y,
# Dim{:azimuth}) 3-D Raster — one layer per direction — so each can be
# resampled independently before being repacked.
function _compute_horizon_native(elevation::Raster; directions::Int)
    canonical = _canonicalize(elevation)
    raw = Geomorphometry.horizon_angle(parent(canonical);
        directions, cellsize = Geomorphometry.cellsize(canonical))
    # Geomorphometry labels slot k as compass (k-1)*(360/n) on the image-matrix
    # convention (row 1 = north, increasing row = south, col 1 = west). Our
    # canonical layout is (X ascending, Y ascending), so row 1 = west and the
    # "N" sweep walks west→east — every slot's physical direction sits 90° CCW
    # of its label. Rotate the azimuth dim by +n/4 (shift indices by -n/4) so
    # slot k once again points at compass (k-1)*(360/n) physically. Without
    # this shift, mid-morning lookups for the "east" horizon return the
    # northern horizon, producing perfect vertical stripes in direct radiation.
    corrected = circshift(raw, (0, 0, -directions ÷ 4))
    az_dim = Dim{:azimuth}(1:directions)
    return Raster(corrected, (dims(canonical, X), dims(canonical, Y), az_dim);
                  crs = crs(canonical))
end

# Resample each azimuth slice of `horizon_native` (X, Y, :azimuth) onto
# `target` with the per-direction mean, then pack into an `SVector{N, °}`
# per pixel.
function _resample_horizon(horizon_native::Raster, target, n::Int)
    target === horizon_native && return _pack_horizon_svector(horizon_native, n)
    xy = dims(target, (X, Y))
    az_dim = Dim{:azimuth}(1:n)
    resampled = Raster(zeros(Float64, length(xy[1]), length(xy[2]), n), (xy..., az_dim))
    @inbounds for k in 1:n
        layer = Rasters.resample(view(horizon_native, Dim{:azimuth}(k));
                                 to = target, method = :average)
        resampled[:, :, k] .= parent(layer)
    end
    return _pack_horizon_svector(resampled, n)
end

# Reinterpret a (X, Y, azimuth) 3-D Raster into a 2-D Raster of
# `SVector{N, °}` so downstream `terrain.horizon_angles[I...]` returns the
# whole horizon vector for one pixel in one indexing step.
function _pack_horizon_svector(src_3d::Raster, n::Int)
    raw = parent(src_3d)
    permuted = permutedims(raw, (3, 1, 2)) .* u"°"
    data = reinterpret(reshape, SVector{n, eltype(permuted)}, permuted)
    xy = dims(src_3d, (X, Y))
    return Raster(collect(data), xy; crs = crs(src_3d))
end

# Isotropic sky-view factor (Dozier & Frew 1990, flat-surface limit):
# V_sky = (1/N) Σ cos²(h_i). `V_sky = 1` for an unobstructed horizon.
@inline _sky_view_from_horizon(horizon_angles) =
    sum(h -> cos(h)^2, horizon_angles) / length(horizon_angles)

# Canonicalise a Raster to `(X, Y)` dim order with both axes ascending.
function _canonicalize(r::Raster)
    r2 = dimnum(r, X) > dimnum(r, Y) ? permutedims(r, (X, Y)) : r
    y_lookup = lookup(r2, Y)
    if length(y_lookup) > 1 && first(y_lookup) > last(y_lookup)
        r2 = reverse(r2; dims = Y)
    end
    return r2
end
