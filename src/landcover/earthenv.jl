# EarthEnv consensus 12-class land cover (fractional).
#
# EarthEnv ships one global GeoTIFF per consensus class — each pixel
# carries a 0-100 percent value for that class — so the loader assembles
# a `RasterStack` keyed by class name. Per-class albedo / roughness
# tables follow the same key set.
#
# ! PLACEHOLDER VALUES — see the warning at the top of `landcover.jl`.

# FIXME: these are guesses
default_landcover_albedo(::Type{<:EarthEnv{<:LandCover}}) = (
    needleleaf_trees = 0.09,
    evergreen_broadleaf_trees = 0.13,
    deciduous_broadleaf_trees = 0.16,
    other_trees = 0.13,
    shrubs = 0.20,
    herbaceous = 0.23,
    cultivated_and_managed = 0.20,
    regularly_flooded = 0.12,
    urban_builtup = 0.15,
    snow_ice = 0.70,
    barren = 0.30,
    open_water = 0.06,
)

# FIXME: these are guesses
default_landcover_roughness(::Type{<:EarthEnv{<:LandCover}}) = (
    needleleaf_trees = 1.0u"m",
    evergreen_broadleaf_trees = 2.0u"m",
    deciduous_broadleaf_trees = 1.5u"m",
    other_trees = 1.2u"m",
    shrubs = 0.10u"m",
    herbaceous = 0.03u"m",
    cultivated_and_managed = 0.05u"m",
    regularly_flooded = 0.05u"m",
    urban_builtup = 1.0u"m",
    snow_ice = 0.001u"m",
    barren = 0.005u"m",
    open_water = 0.0002u"m",
)

"""
    load_landcover(::Type{<:EarthEnv{<:LandCover}}, area::Extent) -> RasterStack

Load every EarthEnv land-cover fractional class for `area`, lazy-read each
tile and crop to `area` before materialising.
"""
function load_landcover(::Type{T}, area::Extent) where {T<:EarthEnv{<:LandCover}}
    # TODO this should just load as a RasterStack directly
    class_names = keys(default_landcover_albedo(T))
    layers = map(class_names) do name
        path = getraster(T, name)
        crop(Raster(path; lazy=true); to=area, touches=true) |> read
    end
    return RasterStack(NamedTuple{class_names}(layers))
end
