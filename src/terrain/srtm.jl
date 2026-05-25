# SRTM 5°×5° tile mosaic.
#
# `Raster(SRTM; extent=…)` (in the dev branch of Rasters) fetches just the
# tiles that intersect `extent` and assembles them into a lazy raster; we
# crop+read locally to materialise the requested window.

# Fractional-degree buffer around the requested area so the resample-to
# template downstream never hits a nodata edge.
const _AREA_BUFFER = 0.05

function _load_dem(::Type{SRTM}, area::Extent)
    buffered = Extents.grow(area, _AREA_BUFFER)
    dem_full = Raster(SRTM; extent = buffered, lazy = true, missingval = Int16(0))
    return read(crop(dem_full; to = buffered, touches = true))
end
