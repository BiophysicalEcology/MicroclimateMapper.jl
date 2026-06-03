# NCEPHourly — 6-hourly NCEP/NCAR Reanalysis bindings.
#
# The NCEP archive stores met variables and radiation/precipitation in different
# dataset categories at different temporal resolutions:
#
#   NCEP{SixHourlyPressure} — air temperature, U/V wind, specific humidity
#                              at pressure levels (lowest level ≈ surface);
#                              1460 steps/year (6-hourly).
#   NCEP{SixHourlySurface}  — surface pressure (:pres_sfc);
#                              1460 steps/year (6-hourly).
#   NCEP{SurfaceGauss}      — downward shortwave/longwave radiation, precip rate;
#                              sub-daily (may be 3-hourly = 2920/year), aggregated
#                              to 6-hourly after loading.
#
# A custom `_load_weather` merges these three sub-stacks before the canonical
# pipeline sees the data. The merged stack has native field names that the
# `WeatherVariable` declarations in `weather_variables(NCEPHourly)` map to the
# `*_6h` staging buffers in `_allocate_weather_buffers(::SixHourlyResolution)`.
#
# The `_DERIVATIONS_6H_TO_1H` chain (see weather.jl) then disaggregates/
# interpolates the staging buffers to hourly canonical output:
#   - Linear interpolation for met variables (T, wind, VP, pressure, LW)
#   - Solar-geometry-aware disaggregation for shortwave radiation
#   - Block summation for rainfall (6h → daily)

"""
    NCEPHourly

Weather source sentinel type for 6-hourly NCEP/NCAR Reanalysis 1 data,
disaggregated to true hourly output via solar-geometry-aware shortwave
disaggregation and linear interpolation for met variables.

Data is loaded from three NCEP dataset categories and merged:
- `NCEP{SixHourlyPressure}` — air temperature, U/V wind, specific humidity
- `NCEP{SixHourlySurface}`  — surface pressure
- `NCEP{SurfaceGauss}`      — shortwave/longwave radiation, precipitation

```julia
model = MicroMapModel(;
    micro_model    = ...,
    dem_source     = SRTM,
    weather_source = NCEPHourly,
)
```
"""
struct NCEPHourly end

temporal_resolution(::Type{NCEPHourly}) = SixHourlyResolution()
weather_derivations(::Type{NCEPHourly}) = _DERIVATIONS_6H_TO_1H

weather_area_buffer(::Type{NCEPHourly}) = 4.0

# 6-hourly NCEP variables — written into the `*_6h` staging buffers.
# Native field names match the keys in the merged stack produced by
# `_load_weather(::Type{NCEPHourly}, ...)` below.
function weather_variables(::Type{NCEPHourly})
    (
        WeatherVariable(:reference_temperature_6h,        :air,     u"K"),
        WeatherVariable(:u_wind_6h,                       :uwnd,    u"m/s"),
        WeatherVariable(:v_wind_6h,                       :vwnd,    u"m/s"),
        WeatherVariable(:specific_humidity_6h,            :shum,    1),
        WeatherVariable(:surface_pressure_6h,             :pres_sfc, u"Pa"),
        WeatherVariable(:downward_shortwave_radiation_6h, :dswrf,   u"W/m^2"),
        WeatherVariable(:longwave_radiation_6h,           :dlwrf,   u"W/m^2"),
        # prate is kg/m²/s; multiply by 6 h × 3600 s/h → kg/m² per 6 h block.
        WeatherVariable(:rainfall_6h,                     :prate,   u"kg/m^2",
                        raw -> raw * 21600.0),
    )
end

# Custom loader: merges three NCEP sub-stacks and normalises radiation/precip
# from whatever sub-daily resolution they arrive at down to 6-hourly.
function _load_weather(::Type{NCEPHourly}, area::Extent, years)
    # 6-hourly met from the lowest pressure level (≈ surface conditions)
    met = _load_layers(YearlyTimeSeries(), NCEP{SixHourlyPressure},
        (:air, :uwnd, :vwnd, :shum), area, years)

    # Surface pressure at 6-hourly resolution
    sfc = _load_layers(YearlyTimeSeries(), NCEP{SixHourlySurface},
        (:pres_sfc,), area, years)

    # Radiation and precipitation from SurfaceGauss — may arrive at 3-hourly
    # (2920/year) and must be aggregated to 6-hourly (1460/year).
    rad = _load_layers(YearlyTimeSeries(), NCEP{SurfaceGauss},
        (:dswrf, :dlwrf, :prate), area, years)

    n6h = 1460 * length(years)
    rad_norm = map((:dswrf, :dlwrf, :prate)) do name
        layer = getproperty(rad, name)
        n = size(layer, Ti)
        n == n6h && return layer
        n % n6h == 0 || error(
            "NCEPHourly: $name has Ti=$n, not a multiple of expected 6h count $n6h")
        return _aggregate_ti_to_daily(layer, n ÷ n6h)
    end

    merged = RasterStack(merge(
        NamedTuple(met),
        NamedTuple(sfc),
        NamedTuple{(:dswrf, :dlwrf, :prate)}(rad_norm),
    ))
    return Rasters.replace_missing(merged, NaN)
end
