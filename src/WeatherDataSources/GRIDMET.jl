@doc """
    get_weather(GRIDMET, point; tstart, tend, kwargs...) -> NamedTuple

Fetch gridMET daily meteorological data for a point location and return a NamedTuple
of Microclimate.jl environment objects plus raw forcing vectors.

gridMET provides daily data at ~4 km resolution for the contiguous United States,
1979 to present. Annual NetCDF files are downloaded and cached via `RasterDataSources.GRIDMET`.

# Arguments
- `point`: GeoInterface-compatible point geometry (e.g. `Point([lon, lat])`); WGS-84
- `tstart`, `tend`: simulation period as `Date` values (full calendar years are downloaded;
  data are filtered to `[tstart, tend]`). Defaults: `tend = tstart`.

# Keyword arguments
- `elevation`: site elevation as a Unitful quantity (e.g. `1800.0u"m"`). Used for
  atmospheric pressure and, together with `grid_elevation`, for lapse-rate correction
  and RH adjustment. Defaults to `nothing` (0 m assumed).
- `grid_elevation`: gridMET grid-cell elevation at `(lon, lat)`. Defaults to `elevation`
  (no lapse correction) if not provided. Provide a WorldClim or SRTM-derived value for
  accurate adjustment (equivalent to the reference used by `micro_usa.R`).
- `albedo`: surface albedo used for the clear-sky solar baseline in cloud-cover estimation
  (default: `0.15`).
- `lapse_rate_type`: temperature lapse rate formulation (default: `EnvironmentalLapseRate()`).
- `vapour_pressure_method`: FluidProperties vapour-pressure equation for RH correction
  (default: `GoffGratch()`).

# Returns
`NamedTuple` with fields:
- `environment_minmax::DailyMinMaxEnvironment` — per-day min/max for `simulate_microclimate`
- `environment_daily::DailyTimeseries`
- `environment_hourly::HourlyTimeseries` — all optional fields set to `nothing`
- `latitude` — site latitude with `u"°"` units
- `days` — day-of-year integers for the requested period
- `tminn`, `tmaxx` — lapse-corrected daily min/max temperature (Unitful K)
- `rhminn`, `rhmaxx` — vapour-pressure-corrected relative humidity (0–1 fractions)
- `ccmax`, `ccmin` — cloud-cover fractions derived from solar radiation (0–1)
- `wind_2m` — wind speed corrected to 2 m height (Unitful m/s)
- `rainfall` — daily precipitation (Unitful kg/m²)
- `tannul` — annual-mean temperature repeated daily (Unitful K)

# Derivation notes
The corrections mirror those in `NicheMapR::micro_usa.R`:
- Temperature: lapse rate applied to tmmx/tmmn (K → Unitful K at site elevation)
- RH: actual vapour pressure conserved; rmin (occurs at Tmax) corrected to site Tmax,
  rmax (occurs at Tmin) corrected to site Tmin
- Wind: power-law height correction from 10 m to 2 m (exponent 0.15); floor 0.1 m/s
- Cloud: `(1 − srad/clearsky)×100`; CCMAX = cloud×2 capped at 1.0, CCMIN = cloud×0.5

# Example
```julia
using BiophysicalGrids, Dates
gm = get_weather(GRIDMET, Point([-106.8, 39.9]);
                 tstart = Date(2010), tend = Date(2012),
                 elevation = 1800.0u"m", grid_elevation = 1750.0u"m")
gm.tminn   # Unitful K daily min temperatures
gm.rhminn  # 0–1 RH fractions
```
""" GRIDMET

function get_weather(
    ::Type{GRIDMET},
    point;
    tstart::Date,
    tend::Date = tstart,
    elevation = nothing,
    grid_elevation = nothing,
    albedo::Float64 = 0.15,
    lapse_rate_type::LapseRate = EnvironmentalLapseRate(),
    vapour_pressure_method = GoffGratch(),
)
    lon, lat = _lonlat(point)

    # ── 1. Download/cache annual files and extract point time series ─────────
    years      = year(tstart):year(tend)
    gm_layers  = (:tmmx, :tmmn, :pr, :rmax, :rmin, :srad, :vs)
    raw        = Dict(l => Float64[] for l in gm_layers)

    for yr in years
        println("  Fetching gridMET $yr...")
        paths = getraster(GRIDMET, gm_layers; date = Date(yr))
        for lyr in gm_layers
            r  = Raster(paths[lyr]; lazy = true)
            ts = r[X(Near(lon)), Y(Near(lat))]
            append!(raw[lyr], Float64.(collect(ts)))
        end
    end

    # ── 2. Date mask (full calendar years downloaded; slice to [tstart, tend]) ─
    all_dates = vcat([collect(Date(yr, 1, 1):Day(1):Date(yr, 12, 31)) for yr in years]...)
    mask  = tstart .<= all_dates .<= tend
    doys  = dayofyear.(all_dates[mask])
    ndays = sum(mask)

    # ── 3. Temperature: lapse-rate correction (stay in Unitful K) ───────────
    site_elev = isnothing(elevation)      ? 0.0u"m" : elevation
    ref_elev  = isnothing(grid_elevation) ? site_elev : grid_elevation
    Δz        = site_elev - ref_elev

    tmaxx_K_grid = raw[:tmmx][mask] .* u"K"
    tminn_K_grid = raw[:tmmn][mask] .* u"K"
    tmaxx_K      = lapse_adjust_temperature(tmaxx_K_grid, Δz, lapse_rate_type)
    tminn_K      = lapse_adjust_temperature(tminn_K_grid, Δz, lapse_rate_type)

    # ── 4. RH: conserve actual vapour pressure ────────────────────────────────
    # rmin (minimum daily RH, occurs at Tmax) corrected from grid Tmax → site Tmax
    # rmax (maximum daily RH, occurs at Tmin) corrected from grid Tmin → site Tmin
    rmin_frac = raw[:rmin][mask] ./ 100.0
    rmax_frac = raw[:rmax][mask] ./ 100.0
    rhminn = clamp.(rh_at_temperature(rmin_frac, tmaxx_K_grid, tmaxx_K, vapour_pressure_method), 0.0001, 1.0)
    rhmaxx = clamp.(rh_at_temperature(rmax_frac, tminn_K_grid, tminn_K, vapour_pressure_method), 0.0001, 1.0)

    # ── 5. Wind: 10 m → 2 m height correction (power law, exponent 0.15) ────
    vs_2m = max.(raw[:vs][mask], 0.1) .* (2.0 / 10.0)^0.15 .* u"m/s"

    # ── 6. Cloud cover from shortwave radiation ───────────────────────────────
    P_atm = atmospheric_pressure(site_elev)
    solar_terrain_flat = SolarTerrain(;
        elevation            = site_elev,
        slope                = 0.0u"°",
        aspect               = 0.0u"°",
        horizon_angles       = fill(0.0u"°", 24),
        albedo               = albedo,
        atmospheric_pressure = P_atm,
        latitude             = lat * u"°",
        longitude            = lon * u"°",
    )
    cloud = cloud_from_srad(raw[:srad][mask] .* u"W/m^2", solar_terrain_flat, doys)
    ccmax = clamp.(cloud .* 2.0, 0.0, 1.0)
    ccmin = clamp.(cloud .* 0.5, 0.0, 1.0)

    # ── 7. Annual mean temperature, rainfall ─────────────────────────────────
    tannul   = fill(mean((tmaxx_K .+ tminn_K) ./ 2), ndays)
    rainfall = raw[:pr][mask] .* u"kg/m^2"

    # ── 8. Standard Microclimate.jl environment objects ──────────────────────
    environment_minmax = DailyMinMaxEnvironment(;
        reference_temperature_min = tminn_K,
        reference_temperature_max = tmaxx_K,
        reference_wind_min        = vs_2m .* 0.1,
        reference_wind_max        = vs_2m,
        reference_humidity_min    = rhminn,
        reference_humidity_max    = rhmaxx,
        cloud_min                 = ccmin,
        cloud_max                 = ccmax,
        minima_times              = fill(0.0, ndays),
        maxima_times              = fill(0.0, ndays),
    )

    environment_daily = DailyTimeseries(;
        shade                 = zeros(ndays),
        soil_wetness          = zeros(ndays),
        surface_emissivity    = fill(0.95, ndays),
        cloud_emissivity      = fill(0.95, ndays),
        rainfall,
        deep_soil_temperature = tannul,
        leaf_area_index       = fill(0.1, ndays),
    )

    environment_hourly = HourlyTimeseries(;
        pressure              = fill(P_atm, ndays * 24),
        reference_temperature = nothing,
        reference_humidity    = nothing,
        reference_wind_speed  = nothing,
        global_radiation      = nothing,
        longwave_radiation    = nothing,
        cloud_cover           = nothing,
        rainfall              = nothing,
        zenith_angle          = nothing,
    )

    return (;
        environment_minmax,
        environment_daily,
        environment_hourly,
        latitude = lat * u"°",
        days     = doys,
        # Raw Unitful forcing vectors for callers that build their own MicroProblem
        tminn    = tminn_K,
        tmaxx    = tmaxx_K,
        rhminn,                # 0–1 fractions
        rhmaxx,
        ccmax,                 # 0–1 fractions
        ccmin,
        wind_2m  = vs_2m,
        rainfall,
        tannul,
    )
end
