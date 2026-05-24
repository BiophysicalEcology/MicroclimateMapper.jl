const _TC_LAYERS = (:tmax, :tmin, :ws, :vap, :vpd, :srad, :ppt, :soil)

# ---------------------------------------------------------------------------
# Per-worker scratch buffers for `assemble_weather!`
# ---------------------------------------------------------------------------

"""
    allocate_weather_buffers(nyears::Int) -> NamedTuple

Per-worker scratch for `assemble_weather!`. Holds every per-pixel array plus
pre-constructed `MonthlyMinMaxEnvironment`/`DailyTimeseries`/`HourlyTimeseries`
referencing those arrays — so an `assemble_weather!` call mutates the buffers
in place and the returned env structs are simply the same objects the buffer
was built with. Allocate once per worker.
"""
function allocate_weather_buffers(nyears::Int)
    nmonths = nyears * 12

    f0(n) = zeros(Float64, n)
    tmax = f0(nmonths); tmin = f0(nmonths); ws = f0(nmonths); vpd = f0(nmonths)
    srad = f0(nmonths); ppt = f0(nmonths); soil = f0(nmonths)
    srad_W = zeros(typeof(0.0u"W/m^2"), nmonths)

    tmax_K = zeros(typeof(0.0u"K"), nmonths)
    tmin_K = zeros(typeof(0.0u"K"), nmonths)
    tmean_K = zeros(typeof(0.0u"K"), nmonths)
    deep_T = zeros(typeof(0.0u"K"), nmonths)

    rh_min = f0(nmonths); rh_max = f0(nmonths)
    rh_min_clamped = f0(nmonths); rh_max_clamped = f0(nmonths)

    ws_max = zeros(typeof(0.0u"m/s"), nmonths)
    ws_min = zeros(typeof(0.0u"m/s"), nmonths)

    cloud = f0(nmonths); cloud_min = f0(nmonths); cloud_max = f0(nmonths)
    soil_moisture_monthly = f0(nmonths)
    ppt_kg = zeros(typeof(0.0u"kg/m^2"), nmonths)
    # `pressure` is constant per pixel (barometric formula on elevation),
    # so the env struct just holds a `Fill` rather than `nmonths*24` copies.
    pressure_placeholder = Fill(atmospheric_pressure(0.0u"m"), nmonths * 24)

    doys = repeat(MONTHLY_BASE_DAYS, nyears)
    days_buf = Float64.(doys)

    # Static (initialised once, never mutated per pixel).
    shade = f0(nmonths); soil_wetness = f0(nmonths)
    surface_emissivity = fill(0.95, nmonths)
    cloud_emissivity   = fill(0.95, nmonths)
    leaf_area_index    = fill(0.1, nmonths)

    environment_minmax = MonthlyMinMaxEnvironment(;
        reference_temperature_min = tmin_K,
        reference_temperature_max = tmax_K,
        reference_wind_min        = ws_min,
        reference_wind_max        = ws_max,
        reference_humidity_min    = rh_min_clamped,
        reference_humidity_max    = rh_max_clamped,
        cloud_min                 = cloud_min,
        cloud_max                 = cloud_max,
        minima_times = (temp = 0, wind = 0, humidity = 1, cloud = 1),
        maxima_times = (temp = 1, wind = 1, humidity = 0, cloud = 0),
    )
    environment_daily = DailyTimeseries(;
        shade, soil_wetness, surface_emissivity, cloud_emissivity,
        rainfall = ppt_kg, deep_soil_temperature = deep_T, leaf_area_index,
    )
    environment_hourly = HourlyTimeseries(;
        pressure              = pressure_placeholder,
        reference_temperature = nothing, reference_humidity   = nothing,
        reference_wind_speed  = nothing, global_radiation     = nothing,
        longwave_radiation    = nothing, cloud_cover          = nothing,
        rainfall              = nothing, zenith_angle         = nothing,
    )

    return (;
        tmax, tmin, ws, vpd, srad, ppt, soil, srad_W,
        tmax_K, tmin_K, tmean_K,
        rh_min, rh_max, rh_min_clamped, rh_max_clamped,
        ws_max, ws_min,
        cloud, cloud_min, cloud_max,
        deep_T, soil_moisture_monthly, ppt_kg,
        doys, days_buf,
        shade, soil_wetness, surface_emissivity, cloud_emissivity, leaf_area_index,
        environment_minmax, environment_daily, environment_hourly,
    )
end

"""
    load_terraclimate(area::Extent, years; scenario=Historical) -> RasterStack

Load every TerraClimate layer needed by the microclimate model over `area`
for `years`, cropped to `area`, concatenated along `Ti` (12 months per year).

Layers: `tmax`, `tmin`, `ws`, `vap`, `vpd`, `srad`, `ppt`, `soil`. The
returned `RasterStack` is the input both to `assemble_weather` for
per-pixel construction of environment structs and to `Rasters.resample`
when aligning weather onto a template grid.
"""
function load_terraclimate(area::Extent, years;
    scenario::Type{<:RasterDataSources.WarmingScenario} = Historical,
)
    layer_rasters = map(_TC_LAYERS) do name
        per_year = map(years) do yr
            path = getraster(TerraClimate{scenario}, name; date = Date(yr))
            read(crop(Raster(path; lazy = true); to = area, touches = true))
        end
        cat(per_year...; dims = Ti)
    end
    return RasterStack(NamedTuple{_TC_LAYERS}(layer_rasters))
end

"""
    assemble_weather!(buf, pixel, elevation; latitude, longitude,
                      solar_out, solar_buffers, kwargs...) -> NamedTuple

Build `(environment_minmax, environment_daily, environment_hourly, ...)` from
one pixel of a `load_terraclimate` stack. `pixel` is a NamedTuple of 1-D
views along `Ti` (one per TerraClimate layer).

`buf` is per-worker scratch from `allocate_weather_buffers` — every
per-pixel array lives there, and the returned env structs are the ones
pre-built by the buffer constructor with their underlying arrays now
holding this pixel's data.

Per-pixel inputs:
- `elevation` — site elevation `Quantity{m}`
- `latitude`, `longitude` — `Quantity{°}` (used for solar geometry only)
- `grid_elevation` — TerraClimate grid-cell elevation for lapse correction
  (defaults to `elevation` ⇒ no correction)
- `vapour_pressure_method`, `lapse_rate_type` — same semantics as
  `get_weather`.
- `solar_out`, `solar_buffers` — per-worker solar scratch (see
  `allocate_output_arrays`, `allocate_buffers`).
"""
function assemble_weather!(scratch, weather, site, I::CartesianIndex{2};
    grid_elevation = site.elevation,
    vapour_pressure_method = GoffGratch(),
    lapse_rate_type::LapseRate = EnvironmentalLapseRate(),
)
    buf = scratch.weather
    @inbounds for k in eachindex(buf.tmax)
        buf.tmax[k] = Float64(weather.tmax[I, k])
        buf.tmin[k] = Float64(weather.tmin[I, k])
        buf.ws[k]   = Float64(weather.ws[I, k])
        buf.vpd[k]  = Float64(weather.vpd[I, k])
        buf.srad[k] = Float64(weather.srad[I, k])
        buf.ppt[k]  = Float64(weather.ppt[I, k])
        buf.soil[k] = Float64(weather.soil[I, k])
    end
    nmonths = length(buf.tmax)
    nyears  = nmonths ÷ 12
    P_atm   = atmospheric_pressure(site.elevation)

    # Temperatures in K; optional lapse correction.
    Δz = site.elevation - grid_elevation
    if iszero(Δz)
        @. buf.tmax_K = u"K"(buf.tmax * u"°C")
        @. buf.tmin_K = u"K"(buf.tmin * u"°C")
    else
        @. buf.tmax_K = lapse_adjust_temperature(u"K"(buf.tmax * u"°C"),
                                                 Δz, lapse_rate_type)
        @. buf.tmin_K = lapse_adjust_temperature(u"K"(buf.tmin * u"°C"),
                                                 Δz, lapse_rate_type)
    end

    # Match NicheMapR micro_terra.R: RH_max paired with T_min, RH_min with T_max.
    @. buf.tmean_K = (buf.tmax_K + buf.tmin_K) / 2
    _rh_from_vpd_at_tmean!(buf.rh_max, buf.vpd, buf.tmean_K, buf.tmin_K,
                           vapour_pressure_method)
    _rh_from_vpd_at_tmean!(buf.rh_min, buf.vpd, buf.tmean_K, buf.tmax_K,
                           vapour_pressure_method)
    @. buf.rh_min_clamped = clamp(buf.rh_min, 0.0, 1.0)
    @. buf.rh_max_clamped = clamp(buf.rh_max, 0.0, 1.0)

    # Wind: 10 m → 2 m via power-law shear (exponent 0.15).
    shear_factor = (2.0 / 10.0)^0.15
    @. buf.ws_max = buf.ws * shear_factor * u"m/s"
    @. buf.ws_min = buf.ws_max * 0.1

    # Cloud cover from observed srad vs flat-terrain clear-sky.
    # `flat_terrain_template` is constant except for elevation, pressure,
    # latitude, longitude — swap them in via `setproperties` (zero-cost).
    @. buf.srad_W = buf.srad * u"W/m^2"
    flat_terrain = setproperties(scratch.cloud_constants.flat_terrain_template,
        (; site.elevation, atmospheric_pressure = P_atm,
           site.latitude, site.longitude))
    cloud_from_srad!(buf.cloud, buf.srad_W, flat_terrain, buf.doys,
                     scratch.solar.out, scratch.solar.buffers;
                     solar_model = scratch.cloud_constants.solar_model,
                     hours       = scratch.cloud_constants.hours,
                     days_buf    = buf.days_buf)
    @. buf.cloud_min = clamp(buf.cloud * 0.5, 0.0, 1.0)
    @. buf.cloud_max = clamp(buf.cloud * 2.0, 0.0, 1.0)

    # Deep soil temperature, soil moisture, rainfall units.
    _annual_means_monthly!(buf.deep_T, buf.tmean_K, nyears)
    @. buf.soil_moisture_monthly = buf.soil / 1000.0 * (1.0 - 1.3 / 2.56)
    @. buf.ppt_kg = buf.ppt * u"kg/m^2"

    # Pressure is constant for this pixel — store it as a `Fill` and swap it
    # into the pre-built hourly env struct via `setproperties` (zero-cost).
    pressure = Fill(P_atm, length(buf.environment_hourly.pressure))
    environment_hourly = setproperties(buf.environment_hourly, (; pressure))

    return (; environment_minmax = buf.environment_minmax,
              environment_daily  = buf.environment_daily,
              environment_hourly,
              site.latitude, days = buf.doys,
              soil_moisture_monthly = buf.soil_moisture_monthly)
end

function assemble_weather(pixel, elevation;
    latitude, longitude,
    grid_elevation = elevation,
    vapour_pressure_method = GoffGratch(),
    lapse_rate_type::LapseRate = EnvironmentalLapseRate(),
    solar_out = nothing,
    solar_buffers = nothing,
)
    tmax = Float64.(collect(pixel.tmax))
    tmin = Float64.(collect(pixel.tmin))
    ws   = Float64.(collect(pixel.ws))
    vpd  = Float64.(collect(pixel.vpd))
    srad = Float64.(collect(pixel.srad))
    ppt  = Float64.(collect(pixel.ppt))
    soil = Float64.(collect(pixel.soil))
    nmonths = length(tmax)
    nyears  = nmonths ÷ 12
    P_atm   = atmospheric_pressure(elevation)

    # Temperatures in K (affine °C cannot be summed); optional lapse correction.
    tmax_K = u"K".(tmax .* u"°C")
    tmin_K = u"K".(tmin .* u"°C")
    Δz = elevation - grid_elevation
    if !iszero(Δz)
        tmax_K = lapse_adjust_temperature(tmax_K, Δz, lapse_rate_type)
        tmin_K = lapse_adjust_temperature(tmin_K, Δz, lapse_rate_type)
    end

    # Match NicheMapR micro_terra.R: RH_max paired with T_min, RH_min with T_max.
    tmean_K = (tmax_K .+ tmin_K) ./ 2
    rh_max  = _rh_from_vpd_at_tmean(vpd, tmean_K, tmin_K, vapour_pressure_method)
    rh_min  = _rh_from_vpd_at_tmean(vpd, tmean_K, tmax_K, vapour_pressure_method)

    # Wind: TerraClimate 10 m → 2 m via power-law shear (exponent 0.15).
    ws_max = ws .* ((2.0 / 10.0)^0.15) .* u"m/s"
    ws_min = ws_max .* 0.1

    # Cloud cover from observed srad vs flat-terrain clear-sky.
    flat_terrain = SolarTerrain(;
        elevation, slope = 0.0u"°", aspect = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24),
        albedo = 0.15, atmospheric_pressure = P_atm,
        latitude, longitude,
    )
    doys   = repeat(MONTHLY_BASE_DAYS, nyears)
    cloud  = if isnothing(solar_out) || isnothing(solar_buffers)
        cloud_from_srad(srad .* u"W/m^2", flat_terrain, doys)
    else
        cloud_from_srad!(srad .* u"W/m^2", flat_terrain, doys, solar_out, solar_buffers)
    end

    # Deep soil temperature: annual mean broadcast to monthly.
    deep_T = _annual_means_monthly(tmean_K, nyears)

    # Soil moisture: TerraClimate column water (mm) → m³/m³ via porosity.
    soil_moisture_monthly = soil ./ 1000.0 .* (1.0 - 1.3 / 2.56)

    environment_minmax = MonthlyMinMaxEnvironment(;
        reference_temperature_min = tmin_K,
        reference_temperature_max = tmax_K,
        reference_wind_min        = ws_min,
        reference_wind_max        = ws_max,
        reference_humidity_min    = clamp.(rh_min, 0.0, 1.0),
        reference_humidity_max    = clamp.(rh_max, 0.0, 1.0),
        cloud_min                 = clamp.(Float64.(cloud) .* 0.5, 0.0, 1.0),
        cloud_max                 = clamp.(Float64.(cloud) .* 2.0, 0.0, 1.0),
        minima_times = (temp = 0, wind = 0, humidity = 1, cloud = 1),
        maxima_times = (temp = 1, wind = 1, humidity = 0, cloud = 0),
    )
    environment_daily = DailyTimeseries(;
        shade              = fill(0.0, nmonths),
        soil_wetness       = fill(0.0, nmonths),
        surface_emissivity = fill(0.95, nmonths),
        cloud_emissivity   = fill(0.95, nmonths),
        rainfall           = ppt .* u"kg/m^2",
        deep_soil_temperature = deep_T,
        leaf_area_index    = fill(0.1, nmonths),
    )
    environment_hourly = HourlyTimeseries(;
        pressure              = fill(P_atm, nmonths * 24),
        reference_temperature = nothing, reference_humidity   = nothing,
        reference_wind_speed  = nothing, global_radiation     = nothing,
        longwave_radiation    = nothing, cloud_cover          = nothing,
        rainfall              = nothing, zenith_angle         = nothing,
    )
    return (; environment_minmax, environment_daily, environment_hourly,
              latitude, days = doys, soil_moisture_monthly)
end

"""
    get_weather(TerraClimate, lon, lat; ystart, yfinish, ...)

Thin per-point wrapper: loads the TerraClimate stack restricted to the
nearest grid cell, then calls `assemble_weather`.
"""
function get_weather(::Type{<:TerraClimate}, lon::Real, lat::Real;
    ystart::Int, yfinish::Int = ystart,
    scenario::Type{<:RasterDataSources.WarmingScenario} = Historical,
    elevation = nothing,
    grid_elevation = nothing,
    lapse_rate_type::LapseRate = EnvironmentalLapseRate(),
    vapour_pressure_method = GoffGratch(),
)
    # Single-pixel area centred on the request point (TerraClimate cell ~4 km).
    pad = 0.05
    area = Extent(X = (lon - pad, lon + pad), Y = (lat - pad, lat + pad))
    stack = load_terraclimate(area, ystart:yfinish; scenario)
    pixel = stack[X = Near(lon), Y = Near(lat)]
    elev  = isnothing(elevation) ? 0.0u"m" : elevation
    gelev = isnothing(grid_elevation) ? elev : grid_elevation
    return assemble_weather(pixel, elev;
        latitude = lat * u"°", longitude = lon * u"°",
        grid_elevation = gelev, vapour_pressure_method, lapse_rate_type,
    )
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""
    _rh_from_vpd_at_tmean!(out, vpd_kPa, Tmean_vec, Tref_vec, method)

In-place version of `_rh_from_vpd_at_tmean`; writes the result into `out`.
"""
function _rh_from_vpd_at_tmean!(
    out::AbstractVector,
    vpd_kPa::AbstractVector,
    Tmean_vec::AbstractVector,
    Tref_vec::AbstractVector,
    method,
)
    @inbounds for k in eachindex(out, vpd_kPa, Tmean_vec, Tref_vec)
        vpd, Tmean, Tref = vpd_kPa[k], Tmean_vec[k], Tref_vec[k]
        e_s_mean_kPa = ustrip(u"hPa", vapour_pressure(method, Tmean)) / 10.0
        actual_e_kPa = e_s_mean_kPa - vpd
        if actual_e_kPa <= 0.0
            out[k] = 0.0
            continue
        end
        e_s_ref_kPa = ustrip(u"hPa", vapour_pressure(method, Tref)) / 10.0
        out[k] = e_s_ref_kPa <= 0.0 ? 1.0 :
                 clamp(actual_e_kPa / e_s_ref_kPa, 0.0, 1.0)
    end
    return out
end

"""
    _rh_from_vpd_at_tmean(vpd_kPa, Tmean_vec, Tref_vec, method) → Vector{Float64}

Derive relative humidity matching the NicheMapR micro_terra.R approach:
  actual_VP = e_sat(Tmean) − VPD
  RH        = actual_VP / e_sat(Tref)

Pass `Tref = Tmin` for RH_max and `Tref = Tmax` for RH_min.
FluidProperties `vapour_pressure` returns hPa; converted to kPa here.
"""
_rh_from_vpd_at_tmean(vpd_kPa, Tmean_vec, Tref_vec, method) =
    _rh_from_vpd_at_tmean!(Vector{Float64}(undef, length(vpd_kPa)),
                           vpd_kPa, Tmean_vec, Tref_vec, method)

"""
    _annual_means_monthly!(out, tmid_all, nyears)

In-place: each month gets that year's mean temperature.
"""
function _annual_means_monthly!(out::AbstractVector, tmid_all::AbstractVector, nyears::Int)
    @inbounds for y in 1:nyears
        s = zero(eltype(tmid_all))
        for m in 1:12
            s += tmid_all[(y - 1) * 12 + m]
        end
        ann_mean = s / 12
        for m in 1:12
            out[(y - 1) * 12 + m] = ann_mean
        end
    end
    return out
end

"""
Broadcast annual mean temperatures to monthly resolution.
For each of `nyears` years (12 months each), every month in the year gets that
year's annual mean temperature. Used for `deep_soil_temperature`.
"""
_annual_means_monthly(tmid_all::AbstractVector, nyears::Int) =
    _annual_means_monthly!(similar(tmid_all), tmid_all, nyears)
