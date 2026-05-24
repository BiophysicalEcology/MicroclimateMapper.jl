"""
    cloud_from_srad(srad, solar_terrain, doy; solar_model, hours)

Estimate fractional cloud cover (0–1) from observed mean daily solar radiation `srad`
(W/m²) by comparing with clear-sky daily mean radiation from `SolarRadiation.solar_radiation`.

`srad` is the observed mean shortwave radiation over 24 hours (including night = 0),
as provided by gridded datasets such as TerraClimate. The clear-sky baseline is computed
using `solar_model` and `solar_terrain` at day-of-year `doy`.

Returns a value clamped to [0, 1]. Returns 1.0 (full cloud cover) if clear-sky radiation
is zero (polar night or numerical issue).

# Example
```julia
cloud = cloud_from_srad(150.0u"W/m^2", solar_terrain, 105)
```
"""
function cloud_from_srad(
    srad,
    solar_terrain::SolarTerrain,
    doy::Int;
    solar_model::SolarProblem = SolarProblem(),
    hours = collect(0.0:1.0:23.0),
)
    sol = solar_radiation(solar_model; solar_terrain, days = [doy], hours)
    clearsky_mean = sum(sol.global_horizontal) / length(hours)
    clearsky_mean <= zero(clearsky_mean) && return 1.0
    cloud = 1.0 - ustrip(u"W/m^2", srad) / ustrip(u"W/m^2", clearsky_mean)
    return clamp(cloud, 0.0, 1.0)
end

"""
    cloud_from_srad(srad_vec, solar_terrain, doys; kwargs...)

Vector version: estimate cloud cover for multiple day-of-year values.
`srad_vec` and `doys` must be the same length.
"""
function cloud_from_srad(
    srad_vec::AbstractVector,
    solar_terrain::SolarTerrain,
    doys::AbstractVector{Int};
    solar_model::SolarProblem = SolarProblem(),
    hours = collect(0.0:1.0:23.0),
)
    # One batched solar_radiation call for all days; allocates internally.
    out = solar_radiation(solar_model; solar_terrain, days = Float64.(doys), hours)
    return _cloud_from_solar_output(srad_vec, out.global_horizontal, length(hours))
end

"""
    cloud_from_srad!(out, srad_vec, solar_terrain, doys, solar_out, solar_buffers; ...)

Fully buffer-reusing version. Caller supplies:
- `out` — cloud-fraction output vector (length `length(doys)`)
- `solar_out` — from `allocate_output_arrays`, sized for `length(hours) * length(doys)`
- `solar_buffers` — from `allocate_buffers`
- `days_buf` (kwarg) — preallocated `Vector{Float64}` for the days argument
  (length `length(doys)`); avoids allocating `Float64.(doys)` per call.
"""
function cloud_from_srad!(
    out::AbstractVector,
    srad_vec::AbstractVector,
    solar_terrain::SolarTerrain,
    doys::AbstractVector{Int},
    solar_out,
    solar_buffers;
    solar_model::SolarProblem = SolarProblem(),
    hours = collect(0.0:1.0:23.0),
    days_buf::AbstractVector{Float64} = Float64.(doys),
)
    @inbounds for k in eachindex(doys, days_buf)
        days_buf[k] = Float64(doys[k])
    end
    solar_radiation!(solar_out, solar_buffers, solar_model;
        solar_terrain, days = days_buf, hours)
    _cloud_from_solar_output!(out, srad_vec, solar_out.global_horizontal, length(hours))
    return out
end

# Allocating wrapper kept for one-off callers.
function cloud_from_srad!(
    srad_vec::AbstractVector,
    solar_terrain::SolarTerrain,
    doys::AbstractVector{Int},
    solar_out,
    solar_buffers;
    kwargs...,
)
    out = Vector{Float64}(undef, length(doys))
    cloud_from_srad!(out, srad_vec, solar_terrain, doys, solar_out, solar_buffers;
        kwargs...)
end

# `global_horizontal` is laid out as ndays blocks of nhours; mean over each
# block gives the per-day clear-sky baseline used to back cloud cover out.
function _cloud_from_solar_output!(out, srad_vec, global_horizontal, nhours)
    ndays = length(srad_vec)
    gh = reshape(global_horizontal, nhours, ndays)
    @inbounds for i in 1:ndays
        clearsky_mean = zero(eltype(global_horizontal))
        for j in 1:nhours
            clearsky_mean += gh[j, i]
        end
        clearsky_mean /= nhours
        if clearsky_mean <= zero(clearsky_mean)
            out[i] = 1.0
        else
            c = 1.0 - ustrip(u"W/m^2", srad_vec[i]) / ustrip(u"W/m^2", clearsky_mean)
            out[i] = clamp(c, 0.0, 1.0)
        end
    end
    return out
end
