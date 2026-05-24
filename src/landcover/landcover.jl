# landcover.jl — load fractional land-cover data and derive per-pixel
# `surface_albedo` and `roughness_height` from it.
#
# Per-class property tables (albedo, roughness length) are **attached to
# the dataset type** via dispatch. EarthEnv's 12-class consensus taxonomy
# is one specific classification; ESA WorldCover, MODIS MCD12, etc. have
# different class sets and would each ship their own table next to their
# own load dispatch.
#
# ! PLACEHOLDER VALUES — NOT YET VERIFIED.
# The numbers in `landcover_albedo_table(::EarthEnv{LandCover})` and
# `landcover_roughness_table(::EarthEnv{LandCover})` below are
# ball-park estimates from general biophysical-modelling knowledge.
# They are NOT pulled from a specific published source. Before relying
# on them for science, replace each value with a citation-backed number
# (broadband shortwave albedo per land-cover class; aerodynamic
# roughness length `z₀` per land-cover class). Otherwise the surface
# parameterisation is junk.

"""
    default_landcover_albedo(::Type{<:LandcoverDataset}) -> NamedTuple

Default per-class broadband surface albedo for the given land-cover
dataset type. Keys must match the layer keys returned by `load_landcover`
for the same dataset. Users can override by passing their own
class→albedo NamedTuple to `microclimate_grid`'s `surface_albedo` kwarg.
"""
function default_landcover_albedo end

"""
    default_landcover_roughness(::Type{<:LandcoverDataset}) -> NamedTuple

Default per-class aerodynamic roughness length `z₀` (length units) for
the given land-cover dataset type. Users can override by passing their
own class→roughness NamedTuple to `microclimate_grid`'s
`roughness_height` kwarg.
"""
function default_landcover_roughness end

# ---------------------------------------------------------------------------
# EarthEnv consensus 12-class land cover
# ---------------------------------------------------------------------------

const _EARTHENV_LANDCOVER_KEYS = (
    :needleleaf_trees, :evergreen_broadleaf_trees, :deciduous_broadleaf_trees,
    :other_trees, :shrubs, :herbaceous, :cultivated_and_managed,
    :regularly_flooded, :urban_builtup, :snow_ice, :barren, :open_water,
)

default_landcover_albedo(::Type{<:EarthEnv{<:LandCover}}) = (
    needleleaf_trees          = 0.09,
    evergreen_broadleaf_trees = 0.13,
    deciduous_broadleaf_trees = 0.16,
    other_trees               = 0.13,
    shrubs                    = 0.20,
    herbaceous                = 0.23,
    cultivated_and_managed    = 0.20,
    regularly_flooded         = 0.12,
    urban_builtup             = 0.15,
    snow_ice                  = 0.70,
    barren                    = 0.30,
    open_water                = 0.06,
)

default_landcover_roughness(::Type{<:EarthEnv{<:LandCover}}) = (
    needleleaf_trees          = 1.0u"m",
    evergreen_broadleaf_trees = 2.0u"m",
    deciduous_broadleaf_trees = 1.5u"m",
    other_trees               = 1.2u"m",
    shrubs                    = 0.10u"m",
    herbaceous                = 0.03u"m",
    cultivated_and_managed    = 0.05u"m",
    regularly_flooded         = 0.05u"m",
    urban_builtup             = 1.0u"m",
    snow_ice                  = 0.001u"m",
    barren                    = 0.005u"m",
    open_water                = 0.0002u"m",
)

"""
    load_landcover(::Type{<:EarthEnv{<:LandCover}}, area::Extent) -> RasterStack

Load every EarthEnv land-cover fractional class for `area`, lazy-read each
tile and crop to `area` before materialising. Layer names match
`_EARTHENV_LANDCOVER_KEYS`.
"""
function load_landcover(::Type{T}, area::Extent) where {T<:EarthEnv{<:LandCover}}
    layers = map(_EARTHENV_LANDCOVER_KEYS) do name
        path = getraster(T, name)
        read(crop(Raster(path; lazy = true); to = area, touches = true))
    end
    return RasterStack(NamedTuple{_EARTHENV_LANDCOVER_KEYS}(layers))
end

# ---------------------------------------------------------------------------
# Per-pixel weighted property derivation
# ---------------------------------------------------------------------------

"""
    landcover_weighted(stack, weights) -> Array

Per-pixel fraction-weighted sum of the layers in `stack`, with `weights[k]`
giving the value to associate with class `k`. Layer fractions are 0–100
(consensus percent), so the result is divided by the per-pixel total to
renormalise (some pixels don't sum to exactly 100).
"""
function landcover_weighted(stack::RasterStack, weights)
    layer_keys = keys(weights)
    layer_arrays = map(name -> parent(stack[name]), layer_keys)
    layer_weights = map(name -> weights[name], layer_keys)
    nx, ny = size(first(layer_arrays))
    out = similar(first(layer_arrays), typeof(first(layer_weights)), nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        total = zero(eltype(first(layer_arrays)))
        acc = zero(typeof(first(layer_weights)))
        ntuple(length(layer_keys)) do k
            f = layer_arrays[k][i, j]
            total += f
            acc += f * layer_weights[k]
            nothing
        end
        out[i, j] = total > 0 ? acc / total : layer_weights[1]
    end
    return out
end

# ---------------------------------------------------------------------------
# Surface property resolution
# ---------------------------------------------------------------------------
#
# `surface_albedo` and `roughness_height` kwargs of `microclimate_grid`
# accept these forms (`prop` below):
#   - `nothing` → use the dataset's default class→property NamedTuple
#     (looked up via `default_fn(landcover_source)`); requires
#     `landcover_source` to be a land-cover dataset type.
#   - a `NamedTuple` keyed by class names → use as the lookup directly
#     (overrides the dataset default).
#   - a scalar → broadcast to every pixel; ignores `landcover_source`.
#   - a `Raster` → resampled to the weather template; ignores
#     `landcover_source`.

# Scalar / unrecognised value: broadcast as a `Fill`.
_resolve_surface_property(value, landcover_source, template, area, default_fn) =
    Fill(value, size(template, X), size(template, Y))

# User-supplied `Raster`: resample to weather template.
_resolve_surface_property(r::Raster, _, template, _, _) =
    parent(Rasters.resample(r; to = template))

# `nothing` → dispatch the default class→property NamedTuple for the
# landcover source, then load + weight.
_resolve_surface_property(::Nothing, landcover_source, template, area, default_fn) =
    _resolve_surface_property(default_fn(landcover_source), landcover_source,
                              template, area, default_fn)

_resolve_surface_property(::Nothing, ::Nothing, _, _, default_fn) =
    error("`$(default_fn)`: no value supplied. Pass a `landcover_source` " *
          "(e.g. `EarthEnv{LandCover}`), or a scalar / Raster / NamedTuple " *
          "for the corresponding surface kwarg.")

# A class→property `NamedTuple` (either the default looked up above or a
# user-supplied override) + a `landcover_source` dataset type.
function _resolve_surface_property(
    weights::NamedTuple, ::Type{T}, template, area, _,
) where {T<:EarthEnv{<:LandCover}}
    stack = load_landcover(T, area)
    weighted_native = landcover_weighted(stack, weights)
    weighted = Raster(weighted_native, dims(first(values(stack))))
    return parent(Rasters.resample(weighted; to = template))
end