# SRTM 5°×5° tile mosaic.
#
# `Raster(SRTM; extent=…)` (in the dev branch of Rasters) fetches just the
# tiles that intersect `extent` and assembles them into a lazy raster; we
# crop+read locally to materialise the requested window.

# Degree buffer around the requested area. Needs to be wide enough that
# horizon-angle ray sweeps from interior pixels never reach the DEM edge
# (otherwise rays terminate early → artificially low horizon angles → too
# little shading at boundary cells). 0.3° ≈ 33 km is well past the
# distance at which a tall mountain still subtends >1° on the horizon.
const _AREA_BUFFER = 0.3

# Additive degrees expansion via `Extents.buffer`. `Extents.grow` is
# multiplicative (no-op on zero-width extents); a single-point input from
# points mode would otherwise hand SRTM a 1×1 DEM and break the
# native-resolution horizon / cellsize calculations.
function _load_dem(::Type{SRTM}, area::Extent)
    buffered = Extents.buffer(area, (X = _AREA_BUFFER, Y = _AREA_BUFFER))
    dem_full = Raster(SRTM; extent = buffered, lazy = true, missingval = Int16(0))
    return read(crop(dem_full; to = buffered, touches = true))
end
