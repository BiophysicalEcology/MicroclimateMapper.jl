# MODIS land-cover bindings.
#
# Currently supports `MODIS{MCD12Q1}` — MODIS Land Cover Type Yearly L3
# Global 500 m product, IGBP 17-class classification (`:LC_Type1`). The
# scheme is categorical: each pixel carries a single integer code in
# `1:17` for its dominant class.
#
# RasterDataSources fetches MODIS via the ORNL DAAC point-API rather than
# as a global gridded product, so `load_landcover` converts an `area::Extent`
# to the API's centre-point + half-widths `(lat, lon, km_ab, km_lr)`
# parameterisation. The API limits each request to ~100 km half-width
# (~200×200 km window). Larger areas would need tiling — not yet supported.
#
# `landcover_class_codes(::Type{MODIS{MCD12Q1}})` declares the IGBP code
# (1-17) for each class name; the generic `landcover_weighted(raster,
# weights, codes)` looks up the right weight per pixel via that mapping.

# IGBP 17-class codes for MCD12Q1's `:LC_Type1` layer.
# Reference: MODIS MCD12Q1 v6 User Guide, IGBP scheme.
landcover_class_codes(::Type{<:MODIS{MCD12Q1}}) = (
    evergreen_needleleaf_forest = 1,
    evergreen_broadleaf_forest = 2,
    deciduous_needleleaf_forest = 3,
    deciduous_broadleaf_forest = 4,
    mixed_forest = 5,
    closed_shrublands = 6,
    open_shrublands = 7,
    woody_savannas = 8,
    savannas = 9,
    grasslands = 10,
    permanent_wetlands = 11,
    croplands = 12,
    urban_builtup = 13,
    cropland_natural_mosaic = 14,
    snow_ice = 15,
    barren = 16,
    water_bodies = 17,
)

# FIXME: these are guesses (same caveat as EarthEnv).
default_landcover_albedo(::Type{<:MODIS{MCD12Q1}}) = (
    evergreen_needleleaf_forest = 0.09,
    evergreen_broadleaf_forest = 0.13,
    deciduous_needleleaf_forest = 0.13,
    deciduous_broadleaf_forest = 0.16,
    mixed_forest = 0.14,
    closed_shrublands = 0.18,
    open_shrublands = 0.22,
    woody_savannas = 0.17,
    savannas = 0.20,
    grasslands = 0.23,
    permanent_wetlands = 0.12,
    croplands = 0.20,
    urban_builtup = 0.15,
    cropland_natural_mosaic = 0.20,
    snow_ice = 0.70,
    barren = 0.30,
    water_bodies = 0.06,
)

# FIXME: these are guesses.
default_landcover_roughness(::Type{<:MODIS{MCD12Q1}}) = (
    evergreen_needleleaf_forest = 1.0u"m",
    evergreen_broadleaf_forest = 2.0u"m",
    deciduous_needleleaf_forest = 1.0u"m",
    deciduous_broadleaf_forest = 1.5u"m",
    mixed_forest = 1.2u"m",
    closed_shrublands = 0.15u"m",
    open_shrublands = 0.05u"m",
    woody_savannas = 0.50u"m",
    savannas = 0.20u"m",
    grasslands = 0.03u"m",
    permanent_wetlands = 0.05u"m",
    croplands = 0.05u"m",
    urban_builtup = 1.0u"m",
    cropland_natural_mosaic = 0.10u"m",
    snow_ice = 0.001u"m",
    barren = 0.005u"m",
    water_bodies = 0.0002u"m",
)

"""
    load_landcover(::Type{<:MODIS{MCD12Q1}}, area::Extent;
                   date::Date = Date(2020, 1, 1)) -> Raster

Fetch the IGBP-coded MCD12Q1 land-cover raster for the given `area`
(converted to MODIS's centre-point + half-width parameterisation) for
the year of `date`. Returns a single 2-D `Raster` of integer class
codes in `1:17`.
"""
function load_landcover(source::Type{<:MODIS{MCD12Q1}}, area::Extent;
                        date::Date = Date(2020, 1, 1))
    # MODIS's ORNL DAAC subset API takes a centre point and half-widths in km
    # (capped at 100) and returns a sinusoidal-grid window. We request the
    # max half-widths so the returned tile reliably encloses `area` after
    # reprojection — sinusoidal-vs-WGS84 warping plus MODIS's pixel snapping
    # mean a smaller request can leave template corners uncovered.
    centre = ((area.X[1] + area.X[2]) / 2, (area.Y[1] + area.Y[2]) / 2)
    path = getraster(source, :LC_Type1;
        lat = centre[2], lon = centre[1], km_ab = 100, km_lr = 100, date)
    return read(Raster(path; lazy = true))
end
