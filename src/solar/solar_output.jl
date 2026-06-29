# solar_output.jl
#
# SolarOutputLayer: user-facing declaration of which solar radiation outputs
# to write alongside (or instead of) microclimate ODE layers.
#
# At init time, `_make_solar_pairs` converts the user's tuple of
# SolarOutputLayer specs into a typed tuple of (spec, wl_range) NamedTuples,
# where wl_range is either `nothing` (broadband) or `NTuple{2,Int}` (contiguous
# wavelength index range). This typed tuple is stored in the cache and iterated
# with `unrolled_map` in the hot pixel loop — no allocation, no runtime dispatch.
#
# Note: values come from scratch.solar.out populated by assemble_weather!.
# These are horizontal-surface radiation at the correct elevation and
# latitude/longitude. Terrain slope/horizon-angle correction of spectral
# output is planned for a future extension.

# ---------------------------------------------------------------------------
# User-facing type
# ---------------------------------------------------------------------------

"""
    SolarOutputLayer(; name, component=:global_horizontal,
                       minimum_wavelength=nothing, maximum_wavelength=nothing)

Declares a solar radiation output layer for `MicroMapModel.solar_output_layers`.

- `name::Symbol` — layer name in the output `RasterStack`
- `component::Symbol` — radiation component:
    - `:global_horizontal`  — total (direct + diffuse) on horizontal surface
    - `:direct_horizontal`  — direct beam on horizontal surface
    - `:diffuse_horizontal` — diffuse sky radiation on horizontal surface
    - `:rayleigh_horizontal`— Rayleigh-scattered component
    - `:global_terrain`     — broadband total on the slope-facing surface
      (broadband only; terrain correction not available for spectral layers)
- `minimum_wavelength::Union{Nothing,Float64}` — lower waveband bound in nm;
  `nothing` (default) uses the pre-integrated broadband scalar.
- `maximum_wavelength::Union{Nothing,Float64}` — upper waveband bound in nm.

The third type parameter `Broadband` is set automatically from the wavelength
arguments — do not supply it manually.

Predefined constants: `SOLAR_BROADBAND`, `SOLAR_PAR`, `SOLAR_UVB`, `SOLAR_NIR`.
"""
struct SolarOutputLayer{Name, Component, Broadband}
    minimum_wavelength::Union{Nothing, Float64}
    maximum_wavelength::Union{Nothing, Float64}
end

function SolarOutputLayer(;
    name::Symbol,
    component::Symbol = :global_horizontal,
    minimum_wavelength::Union{Nothing, Real} = nothing,
    maximum_wavelength::Union{Nothing, Real} = nothing,
)
    broadband = minimum_wavelength === nothing && maximum_wavelength === nothing
    SolarOutputLayer{name, component, broadband}(
        isnothing(minimum_wavelength) ? nothing : Float64(minimum_wavelength),
        isnothing(maximum_wavelength) ? nothing : Float64(maximum_wavelength),
    )
end

@inline _solar_layer_name(::SolarOutputLayer{N}) where N = N
@inline _solar_component(::SolarOutputLayer{N, C}) where {N, C} = C

# ---------------------------------------------------------------------------
# Predefined waveband constants
# ---------------------------------------------------------------------------

const SOLAR_BROADBAND = SolarOutputLayer(; name = :global_radiation)
const SOLAR_PAR       = SolarOutputLayer(; name = :par,
    minimum_wavelength = 400.0, maximum_wavelength = 700.0)
const SOLAR_UVB       = SolarOutputLayer(; name = :uv_b,
    minimum_wavelength = 290.0, maximum_wavelength = 315.0)
const SOLAR_NIR       = SolarOutputLayer(; name = :nir,
    minimum_wavelength = 700.0, maximum_wavelength = 4000.0)

const _DEFAULT_SOLAR_LAYERS = (
    SolarOutputLayer(; name = :global_horizontal),
    SolarOutputLayer(; name = :direct_horizontal,  component = :direct_horizontal),
    SolarOutputLayer(; name = :diffuse_horizontal, component = :diffuse_horizontal),
    SolarOutputLayer(; name = :global_terrain,     component = :global_terrain),
)

# canonical_unit for solar output layers (:global_radiation is in common.jl).
canonical_unit(::Val{:global_horizontal})   = u"W/m^2"
canonical_unit(::Val{:direct_horizontal})   = u"W/m^2"
canonical_unit(::Val{:diffuse_horizontal})  = u"W/m^2"
canonical_unit(::Val{:global_terrain})      = u"W/m^2"
canonical_unit(::Val{:rayleigh_horizontal}) = u"W/m^2"
canonical_unit(::Val{:par})                 = u"W/m^2"
canonical_unit(::Val{:uv_b})                = u"W/m^2"
canonical_unit(::Val{:nir})                 = u"W/m^2"

# ---------------------------------------------------------------------------
# Component → field name dispatch (compile-time, no Dict lookup)
# ---------------------------------------------------------------------------

@inline _solar_broadband_field(::Val{:global_horizontal})   = :global_horizontal
@inline _solar_broadband_field(::Val{:direct_horizontal})   = :direct_horizontal
@inline _solar_broadband_field(::Val{:diffuse_horizontal})  = :diffuse_horizontal
@inline _solar_broadband_field(::Val{:rayleigh_horizontal}) = :rayleigh_horizontal
@inline _solar_broadband_field(::Val{:global_terrain})      = :global_terrain

@inline _solar_spectra_field(::Val{:global_horizontal})   = :global_spectra
@inline _solar_spectra_field(::Val{:direct_horizontal})   = :direct_spectra
@inline _solar_spectra_field(::Val{:diffuse_horizontal})  = :diffuse_spectra
@inline _solar_spectra_field(::Val{:rayleigh_horizontal}) = :rayleigh_spectra

# ---------------------------------------------------------------------------
# Waveband index range computation (runs once at init, not per pixel)
# ---------------------------------------------------------------------------

# Resolve the effective solar layers: when solar_only=true and the user left
# solar_output_layers empty, default to the four standard broadband outputs.
_effective_solar_layers(model) =
    model.solar_only && isempty(model.solar_output_layers) ?
        _DEFAULT_SOLAR_LAYERS : model.solar_output_layers

# For a broadband spec, the range is nothing (use pre-integrated scalar).
_layer_wl_range(::SolarOutputLayer{N, C, true}, _) where {N, C} = nothing

# For a waveband spec, find the contiguous index range in wavelengths.
function _layer_wl_range(layer::SolarOutputLayer{N, C, false}, wavelengths) where {N, C}
    first_index = findfirst(≥(layer.minimum_wavelength * u"nm"), wavelengths)
    last_index  = findlast(≤(layer.maximum_wavelength * u"nm"), wavelengths)
    (first_index === nothing || last_index === nothing || first_index > last_index) && error(
        "SolarOutputLayer :$N: no wavelengths found in " *
        "[$(layer.minimum_wavelength), $(layer.maximum_wavelength)] nm")
    return (first_index, last_index)
end

# Build the typed (spec, wavelength_range) pair NamedTuple for one layer.
# wavelength_range is Nothing or NTuple{2,Int} — a concrete type per pair, so the
# resulting tuple carries exact per-element types and unrolled_map is type-stable.
_solar_layer_pair(spec, wavelengths) =
    (; spec, wavelength_range = _layer_wl_range(spec, wavelengths))

# Build the full typed tuple of pairs for all solar output layers.
_make_solar_pairs(layers::Tuple, wavelengths) =
    map(l -> _solar_layer_pair(l, wavelengths), layers)

# ---------------------------------------------------------------------------
# Per-pixel write helpers
# ---------------------------------------------------------------------------

# Trapezoidal integration over a contiguous wavelength index range.
# spectra_row: view of solar_out.global_spectra[t, :] etc. (W/nm/m²)
# wavelengths: the solar model's wavelength vector (Unitful, nm)
# Returns W/m² — units flow through Unitful.
function _integrate_waveband(spectra_row, first_index::Int, last_index::Int, wavelengths)
    result = 0.0u"W/m^2"
    for k in first_index:(last_index - 1)
        dλ = wavelengths[k + 1] - wavelengths[k]
        result += (spectra_row[k] + spectra_row[k + 1]) * (dλ * 0.5)
    end
    return result
end

# Broadband: read the pre-integrated scalar directly.
# cloud_factors: per-hour Ångström sunshine fractions (nothing = clear-sky, no scaling).
@inline function _write_solar_slice!(rast, solar_out,
                                     ::SolarOutputLayer{N, C, true}, ::Nothing, _, I,
                                     cloud_factors) where {N, C}
    src = getproperty(solar_out, _solar_broadband_field(Val(C)))
    @inbounds for t in eachindex(src)
        rast[I..., Ti(t)] = cloud_factors === nothing ? src[t] : src[t] * cloud_factors[t]
    end
    return nothing
end

# Waveband: trapezoidal integration over the spectral matrix.
@inline function _write_solar_slice!(rast, solar_out,
                                     ::SolarOutputLayer{N, C, false},
                                     wavelength_range::NTuple{2, Int}, wavelengths, I,
                                     cloud_factors) where {N, C}
    spectra = getproperty(solar_out, _solar_spectra_field(Val(C)))
    first_index, last_index = wavelength_range
    @inbounds for t in axes(spectra, 1)
        val = _integrate_waveband(view(spectra, t, :), first_index, last_index, wavelengths)
        rast[I..., Ti(t)] = cloud_factors === nothing ? val : val * cloud_factors[t]
    end
    return nothing
end

# Write all solar output layers for pixel I.
# solar_pairs is a typed tuple of (; spec, wavelength_range) NamedTuples — iterated
# with unrolled_map so each pair specialises to a distinct concrete type,
# giving compile-time dispatch to the correct _write_solar_slice! method.
# cloud_factors: per-hour Ångström sunshine fractions, or nothing for clear-sky.
@inline function _write_solar_output!(output, solar_out, solar_pairs::Tuple, wavelengths, I;
                                       cloud_factors = nothing)
    unrolled_map(solar_pairs) do pair
        _write_solar_slice!(output[_solar_layer_name(pair.spec)],
                            solar_out, pair.spec, pair.wavelength_range, wavelengths, I,
                            cloud_factors)
        nothing
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Output allocation
# ---------------------------------------------------------------------------

# Allocate one (spatial…, Ti) Raster per solar layer, element type W/m².
# Ti axis matches the hourly microclimate axis built from anchor_dates.
function _allocate_solar_output(terrain, solar_pairs::Tuple, anchor_dates, mask)
    spatial_dims = dims(terrain.elevation)
    ti = Ti(_ti_datetime_axis(anchor_dates, collect(0.0:1.0:23.0)))
    T = typeof(NaN * u"W/m^2")
    return RasterStack(NamedTuple(map(solar_pairs) do pair
        ds = (spatial_dims..., ti)
        _solar_layer_name(pair.spec) => _allocate_layer_array(T, ds, mask)
    end))
end
