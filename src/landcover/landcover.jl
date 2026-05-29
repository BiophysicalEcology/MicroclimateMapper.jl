# landcover.jl — load fractional or categorical land-cover data and derive
# per-pixel `surface_albedo` and `roughness_height` from it.
#
# Per-class property tables (albedo, roughness length) are attached to the
# dataset type via dispatch. EarthEnv's 12-class consensus taxonomy is one
# specific classification; MODIS's MCD12Q1 (IGBP 17-class) is another;
# ESA WorldCover, etc. each ship their own table next to their own loader.
#
# Two physical layouts are supported:
#   * **Fractional**  — a `RasterStack` with one 0-100 layer per class
#                       (EarthEnv consensus, MOD44B vegetation continuous fields).
#   * **Categorical** — a single `Raster` of integer class codes per pixel
#                       (MODIS MCD12Q1). Categorical sources also declare a
#                       `landcover_class_codes` mapping from class name to
#                       integer code.
#
# Dispatch on the loaded landcover's type (`RasterStack` vs `Raster`)
# selects the right weighted-aggregation path inside `_resolve_surface_property`.
#
# ! PLACEHOLDER VALUES — NOT YET VERIFIED.
# The albedo and roughness numbers in this file and in `modis.jl` are
# ball-park estimates from general biophysical-modelling knowledge, not
# from a specific published source. Replace with citation-backed values
# (broadband shortwave albedo per class; aerodynamic roughness length `z₀`
# per class) before using for science.

# ---------------------------------------------------------------------------
# Source-extension points
# ---------------------------------------------------------------------------

"""
    default_landcover_albedo(::Type{<:LandcoverDataset}) -> NamedTuple

Default per-class broadband surface albedo for the given land-cover
dataset type. Keys must match the class names of that dataset (also used
as the keys of `landcover_class_codes` for categorical sources). Users
can override by passing their own class→albedo NamedTuple to
`microclimate_grid`'s `surface_albedo` kwarg.
"""
function default_landcover_albedo end

"""
    default_landcover_roughness(::Type{<:LandcoverDataset}) -> NamedTuple

Default per-class aerodynamic roughness length `z₀` (length units) for
the given land-cover dataset type. Same keying as
`default_landcover_albedo`. Users can override by passing their own
class→roughness NamedTuple to `microclimate_grid`'s `roughness_height`
kwarg.
"""
function default_landcover_roughness end

"""
    load_landcover(::Type{<:LandcoverDataset}, area::Extent)
        -> RasterStack OR Raster

Per-source loader. Fractional sources return a `RasterStack` with one
0-100 layer per class (keyed by class name). Categorical sources return
a single `Raster` of integer class codes — those sources must also
implement `landcover_class_codes`.
"""
function load_landcover end

"""
    landcover_class_codes(::Type{<:LandcoverDataset}) -> NamedTuple

For categorical sources (where `load_landcover` returns a single
integer-coded `Raster`), declare the mapping from class name (the key
used in albedo/roughness NamedTuples) to the integer code that appears
in the loaded raster. Not used by fractional sources.
"""
function landcover_class_codes end

# ---------------------------------------------------------------------------
# Per-pixel weighted property derivation
# ---------------------------------------------------------------------------

"""
    landcover_weighted(stack::RasterStack, weights::NamedTuple, _source) -> Raster

Per-pixel fraction-weighted sum of the layers in `stack`. `weights[name]`
gives the value to associate with class `name`; layer fractions are 0-100
(consensus percent), so the result is renormalised by the per-pixel total
(some pixels don't sum to exactly 100). The `source` argument is ignored —
present only so categorical and fractional dispatch share one signature.
"""
function landcover_weighted(stack::RasterStack{K}, weights::NamedTuple, _source) where K
    # Reorder/subset weights to match the stack's layer order so per-pixel
    # NamedTuples and weight NamedTuples align positionally.
    aligned_weights = NamedTuple{K}(weights)
    fallback = first(aligned_weights)
    return map(stack) do fractions
        total = sum(fractions)
        iszero(total) ? fallback : sum(map(*, fractions, aligned_weights)) / total
    end
end
"""
    landcover_weighted(raster::Raster, weights::NamedTuple, source) -> Raster

Per-pixel lookup for a categorical (integer-coded) landcover raster.
`landcover_class_codes(source)[name]` gives the integer code that represents
class `name` in `raster`; `weights[name]` gives the value to assign to that
class. Pixels whose code isn't in the mapping get the first weight as a
fallback.
"""
function landcover_weighted(raster::Raster, weights::NamedTuple, source)
    codes = landcover_class_codes(source)
    fallback = first(values(weights))
    # Code→weight as a flat array indexed by `code + 1` (codes can be 0).
    max_code = maximum(values(codes))
    weight_by_code = fill(fallback, max_code + 1)
    for name in keys(weights)
        haskey(codes, name) || continue
        weight_by_code[codes[name] + 1] = weights[name]
    end
    return map(raster) do value
        code = Int(value)
        (0 <= code <= max_code) ? weight_by_code[code + 1] : fallback
    end
end

# ---------------------------------------------------------------------------
# Surface property resolution
# ---------------------------------------------------------------------------
#
# The `surface_albedo_source` / `roughness_height_source` fields on
# `MicroMapModel` accept these forms:
#   - `nothing` → use the dataset's default class→property NamedTuple
#     (looked up via `default_fn(landcover_source)`); requires
#     `landcover_source` to be a land-cover dataset type.
#   - a `NamedTuple` keyed by class names → use as the lookup directly
#     (overrides the dataset default).
#   - a scalar → broadcast to every pixel; ignores `landcover_source`.
#   - a `Raster` → resampled to the weather template; ignores
#     `landcover_source`.
#   - a `Type{<:RasterDataSource}` → load a property raster directly from
#     that dataset, resample. No landcover weighting.

# Scalar / unrecognised value: broadcast as a `Fill` wrapped in a Raster
# carrying the template's X/Y dims so downstream `[X(i), Y(j)]` indexing
# stays uniform across all branches.
_resolve_surface_property(value, landcover_source, template, area, default_fn) =
    Raster(Fill(value, size(template, X), size(template, Y)), dims(template, (X, Y)))
# User-supplied `Raster`: resample to weather template.
_resolve_surface_property(r::Raster, _, template, _, _) =
    Rasters.resample(r; to=template)
# `nothing` → dispatch the default class→property NamedTuple for the
# landcover source, then load + weight.
_resolve_surface_property(::Nothing, landcover_source, template, area, default_fn) =
    _resolve_surface_property(default_fn(landcover_source), landcover_source,
                              template, area, default_fn)
_resolve_surface_property(::Nothing, ::Nothing, _, _, default_fn) =
    error("`$(default_fn)`: no value supplied. Pass a `landcover_source` " *
          "(e.g. `EarthEnv{LandCover}`, `MODIS{MCD12Q1}`), a `Type{<:RasterDataSource}` " *
          "for the property itself, or a scalar / Raster / NamedTuple.")
# A class→property `NamedTuple` + any landcover source type. Dispatch is
# uniform: load whatever `load_landcover` returns, hand it to the
# resolution-specific weighted aggregator, wrap as a Raster, resample.
function _resolve_surface_property(
    weights::NamedTuple,
    source::Type{<:RasterDataSources.RasterDataSource},
    template, area, _,
)
    landcover = load_landcover(source, area)
    weighted_native = landcover_weighted(landcover, weights, source)
    return _resample_to_template(weighted_native, commondims(landcover, (X(), Y())), template)
end
# A `Type{<:RasterDataSource}` is the property *source* itself (not
# landcover): load a 2-D Raster from that dataset and resample. Ignores
# `landcover_source`. Each dataset that supports being used this way
# implements `load_surface_property(source, area)`.
function _resolve_surface_property(
    source::Type{<:RasterDataSources.RasterDataSource},
    _landcover_source, template, area, _default_fn,
)
    raster = load_surface_property(source, area)
    return Rasters.resample(raster; to = template)
end

"""
    load_surface_property(::Type{<:RasterDataSource}, area::Extent) -> Raster

Per-source loader for surface properties (broadband albedo, roughness
length, …) used when a `MicroMapModel.surface_albedo_source` or
`.roughness_height_source` is set to a `Type{<:RasterDataSource}` directly
(no landcover weighting). The returned `Raster` must be unit-tagged in the
canonical units (`Float64` for albedo, `u"m"` for roughness length).
"""
function load_surface_property end

# Resample a weighted-property array onto the weather `template`. GDAL
# can't write Unitful `Quantity` arrays, so for unit-tagged inputs (e.g.
# roughness in `u"m"`) we strip the unit, resample, then re-attach.
# Returns a `Raster` so downstream callers can index with dim wrappers.
@inline _resample_to_template(values, spatial_dims, template) =
    Rasters.resample(Raster(values, spatial_dims); to = template)
@inline function _resample_to_template(values::AbstractArray{<:Quantity}, spatial_dims, template)
    value_unit = unit(eltype(values))
    stripped = ustrip.(value_unit, values)
    resampled = Rasters.resample(Raster(stripped, spatial_dims); to = template)
    return rebuild(resampled, parent(resampled) .* value_unit)
end

