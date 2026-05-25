"""
    cloud_from_solar_radiation!(out, solar_out, solar_radiation, solar_terrain,
                                days_of_year, solar_buffers; solar_model, hours)

Estimate fractional cloud cover (0–1) for each day in `days_of_year` from
observed mean daily shortwave radiation `solar_radiation[i]` by comparing
against the clear-sky daily mean computed via `SolarRadiation.solar_radiation!`.

`solar_radiation` is the observed mean over 24 hours (including night = 0),
as provided by gridded datasets such as TerraClimate / WorldClim. The result
is written into `out` in-place and clamped to [0, 1]; days with zero
clear-sky radiation (polar night) get `1.0`.

Caller-supplied buffers (written-to objects come first):
- `out` — length `length(days_of_year)`
- `solar_out`     — from `SolarRadiation.allocate_output_arrays`, sized for
                    `length(hours) * length(days_of_year)`
- `solar_buffers` — from `SolarRadiation.allocate_buffers`
"""
function cloud_from_solar_radiation!(
    out::AbstractVector,
    solar_out,
    solar_radiation::AbstractVector,
    solar_terrain::SolarTerrain,
    days_of_year::AbstractVector{Int},
    solar_buffers;
    solar_model::SolarProblem = SolarProblem(),
    hours = collect(0.0:1.0:23.0),
)
    solar_radiation!(solar_out, solar_buffers, solar_model;
        solar_terrain, days = days_of_year, hours)
    _cloud_from_solar_output!(out, solar_radiation,
        solar_out.global_horizontal, length(hours))
    return out
end

# `global_horizontal` is laid out as ndays blocks of nhours; mean over each
# block gives the per-day clear-sky baseline used to back cloud cover out.
function _cloud_from_solar_output!(out, solar_radiation, global_horizontal, nhours)
    ndays = length(solar_radiation)
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
            c = 1.0 - ustrip(u"W/m^2", solar_radiation[i]) / ustrip(u"W/m^2", clearsky_mean)
            out[i] = clamp(c, 0.0, 1.0)
        end
    end
    return out
end
