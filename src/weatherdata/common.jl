const MONTHLY_BASE_DAYS = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]

"""
    monthly_days(nyears) → Vector{Int}

Day-of-year values for the midpoint of each month, repeated for `nyears` years.
Used to set the `days` field of `MicroModel` for monthly forcing data.
"""
monthly_days(nyears::Int) = repeat(MONTHLY_BASE_DAYS, nyears)

"""
    mean_annual_temperature(tmax_vec, tmin_vec) → Quantity

Compute annual mean temperature from vectors of monthly max and min temperatures.
Used as an estimate of deep soil temperature.
"""
function mean_annual_temperature(tmax_vec::AbstractVector, tmin_vec::AbstractVector)
    return mean((tmax_vec .+ tmin_vec) ./ 2)
end

"""
    extract_point(raster, lon, lat)

Extract a scalar (or time-series) value from `raster` at the nearest grid point
to (`lon`, `lat`). Returns the bare array when the raster has a time dimension.
"""
function extract_point(raster, lon::Real, lat::Real)
    return raster[X(Near(lon)), Y(Near(lat))]
end

# ---------------------------------------------------------------------------
# Generic simulation entry point
# ---------------------------------------------------------------------------

"""
    simulate_microclimate(model::MicroModel, site::Site, weather::NamedTuple;
                          clearsky=false,
                          initial_soil_temperature=nothing,
                          initial_soil_moisture=fill(0.42 * 0.25, length(model.depths)),
                          kwargs...)

Run the microclimate model for a single location.

`model` is a fully-built `MicroModel` carrying the simulation recipe (geometry,
physical-process models, soil properties / hydraulics, snow, solver/strategy
choices in `config`). It is constant across calls — swap `site` or `weather`
to run the same model on different data.

`weather` is a `NamedTuple` returned by `get_weather`, containing:
- `environment_minmax::MonthlyMinMaxEnvironment`
- `environment_daily::DailyTimeseries`
- `environment_hourly::HourlyTimeseries` (or `nothing`)

Extra `kwargs...` are forwarded to `MicroInputs`.

# Keyword arguments
- `clearsky`: override `environment_minmax` cloud cover to zero for all
  timesteps (default: `false`).
- `initial_soil_temperature`: initial soil temperatures. Defaults to `nothing`,
  which lets Microclimate.jl reset T₀ to the day-mean reference air temperature
  (matching NicheMapR's monthly `microdaily=0` behaviour).
- `initial_soil_moisture`: initial volumetric soil moisture fractions, sized to
  `model.depths`.

# Returns
`MicroResult` from `Microclimate.solve`.
"""
function simulate_microclimate(
    model::MicroModel,
    site::Site,
    weather::NamedTuple;
    clearsky::Bool = false,
    initial_soil_temperature = nothing,
    initial_soil_moisture = fill(0.42 * 0.25, length(model.depths)),
    kwargs...,
)
    (; environment_minmax, environment_daily, environment_hourly) = weather

    if clearsky
        n = length(environment_minmax.cloud_min)
        environment_minmax = MonthlyMinMaxEnvironment(;
            reference_temperature_min = environment_minmax.reference_temperature_min,
            reference_temperature_max = environment_minmax.reference_temperature_max,
            reference_wind_min        = environment_minmax.reference_wind_min,
            reference_wind_max        = environment_minmax.reference_wind_max,
            reference_humidity_min    = environment_minmax.reference_humidity_min,
            reference_humidity_max    = environment_minmax.reference_humidity_max,
            cloud_min                 = fill(0.0, n),
            cloud_max                 = fill(0.0, n),
            minima_times              = environment_minmax.minima_times,
            maxima_times              = environment_minmax.maxima_times,
        )
    end

    inputs = MicroInputs(;
        site,
        environment_minmax,
        environment_daily,
        environment_hourly,
        initial_soil_temperature,
        initial_soil_moisture,
        kwargs...,
    )

    return Microclimate.solve(MicroProblem(model, inputs))
end
