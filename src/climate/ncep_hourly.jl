# NCEPHourly — 6-hourly NCEP/NCAR Reanalysis bindings.
#
# The NCEP archive stores met variables and radiation/precipitation in different
# dataset categories at different temporal resolutions:
#
#   NCEP{Pressure,1}    — air temperature, U/V wind, specific humidity
#                          at pressure levels (lowest level ≈ surface).
#   NCEP{Surface,1}     — surface pressure (:pres_sfc).
#   NCEP{SurfaceFlux,1} — downward shortwave/longwave radiation, precip rate.
# All are natively 6-hourly (1460 steps/year).
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

Data is loaded from three NCEP groups and merged:
- `NCEP{Pressure,1}`    — air temperature, U/V wind, specific humidity
- `NCEP{Surface,1}`     — surface pressure
- `NCEP{SurfaceFlux,1}` — shortwave/longwave radiation, precipitation

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

# Custom loader: merges the three 6-hourly NCEP groups into one stack.
function _load_weather(::Type{NCEPHourly}, area::Extent, years)
    # 6-hourly met from the lowest pressure level (≈ surface conditions)
    met = _load_layers(YearlyTimeSeries(), NCEP{Pressure,1},
        (:air, :uwnd, :vwnd, :shum), area, years)

    # Surface pressure
    sfc = _load_layers(YearlyTimeSeries(), NCEP{Surface,1},
        (:pres_sfc,), area, years)

    # Radiation and precipitation from the surface-flux group
    rad = _load_layers(YearlyTimeSeries(), NCEP{SurfaceFlux,1},
        (:dswrf, :dlwrf, :prate), area, years)

    merged = RasterStack(merge(
        NamedTuple(met),
        NamedTuple(sfc),
        NamedTuple(rad),
    ))
    return Rasters.replace_missing(merged, NaN)
end
