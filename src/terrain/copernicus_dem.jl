# Copernicus DEM (GLO-30 / GLO-90) 1°×1° Cloud-Optimized GeoTIFF tiles.
#
# `Raster(CopernicusDEM; extent=…, res=…)` fetches the tiles intersecting
# `extent` from the AWS Open Data registry and mosaics them. Absent oceanic
# tiles are filled with `missingval` (0 m), matching the ocean handling of
# `_load_dem(::Type{SRTM}, …)`; we crop+read locally to materialise the window.
#
# `res` selects the resolution ("30m" GLO-30 by default, or "90m" GLO-90).
# `_AREA_BUFFER` (defined for SRTM) is the shared horizon-ray edge buffer.
function _load_dem(::Type{CopernicusDEM}, area::Extent; res = "30m")
    buffered = Extents.buffer(area, (X = _AREA_BUFFER, Y = _AREA_BUFFER))
    dem_full = Raster(CopernicusDEM; extent = buffered, res, lazy = true, missingval = 0.0f0)
    return read(crop(dem_full; to = buffered, touches = true))
end
