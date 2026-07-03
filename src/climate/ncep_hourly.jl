# NCEP — 6-hourly NCEP/NCAR Reanalysis data, disaggregated to hourly.
#
# Currently hardcoded to Reanalysis 1 (NCEP/NCAR, 1948–present) via
# NCEP{SurfaceFlux, 1}. Reanalysis 2 support (NCEP/DOE, 1979–present) would
# follow the same pattern with NCEP{SurfaceFlux, 2}.
#
# All variables are loaded from NCEP{SurfaceFlux, 1} (the T62 Gaussian grid,
# ~1.875°), which provides 6-hourly 2m temperature, 10m wind, 2m specific
# humidity, surface pressure, radiation, and precipitation in a single
# consistent group. This avoids the much larger pressure-level files
# (NCEP{Pressure, 1}) which contain all 17 standard levels when only the
# surface-adjacent level is needed.
#
# Some archive years store radiation/precip at 3-hourly resolution (2920
# steps/year). These are averaged to 6-hourly here to fit the `*_6h` staging
# buffers. The follow-up SubDailyResolution{N} refactor will instead
# interpolate 3-hourly data directly to hourly, preserving full resolution.
#
# The `_DERIVATIONS_6H_TO_1H` chain (see weather.jl) then disaggregates/
# interpolates the staging buffers to hourly canonical output:
#   - Linear interpolation for met variables (T, wind, VP, pressure, LW)
#   - Solar-geometry-aware disaggregation for shortwave radiation
#   - Block summation for rainfall (6h → daily)

"""
    NCEP

Weather source for 6-hourly NCEP/NCAR Reanalysis data disaggregated to
true hourly output via solar-geometry-aware shortwave disaggregation and
linear interpolation for met variables. Currently uses Reanalysis 1
(NCEP/NCAR, 1948–present).

All variables are loaded from `NCEP{SurfaceFlux, 1}` (T62 Gaussian grid,
~1.875°): 2m air temperature, 10m wind, 2m specific humidity, surface
pressure, shortwave/longwave radiation, and precipitation.

```julia
model = MicroMapModel(;
    micro_model    = ...,
    dem_source     = SRTM,
    weather_source = NCEP,
)
```
"""
struct NCEP end

temporal_resolution(::Type{NCEP}) = SixHourlyResolution()
weather_derivations(::Type{NCEP}) = _DERIVATIONS_6H_TO_1H

weather_area_buffer(::Type{NCEP}) = 4.0

# 6-hourly NCEP variables — written into the `*_6h` staging buffers.
# Native field names are the layer keys in NCEP{SurfaceFlux, 1}.
function weather_variables(::Type{NCEP})
    (
        WeatherVariable(:reference_temperature_6h,        :air_2m,   u"K"),
        WeatherVariable(:u_wind_6h,                       :uwnd_10m, u"m/s"),
        WeatherVariable(:v_wind_6h,                       :vwnd_10m, u"m/s"),
        WeatherVariable(:specific_humidity_6h,            :shum_2m,  1),
        WeatherVariable(:surface_pressure_6h,             :pres,     u"Pa"),
        WeatherVariable(:downward_shortwave_radiation_6h, :dswrf,    u"W/m^2"),
        WeatherVariable(:longwave_radiation_6h,           :dlwrf,    u"W/m^2"),
        # prate is kg/m²/s; multiply by 6 h × 3600 s/h → kg/m² per 6 h block.
        WeatherVariable(:rainfall_6h,                     :prate,    u"kg/m^2",
                        raw -> raw * 21600.0),
    )
end

# Load all variables from a single NCEP{SurfaceFlux, 1} stack, then normalise
# any layers that arrive at 3-hourly resolution down to 6-hourly.
function _load_weather(::Type{NCEP}, area::Extent, years)
    stack = _load_layers(YearlyTimeSeries(), NCEP{SurfaceFlux, 1},
        (:air_2m, :shum_2m, :uwnd_10m, :vwnd_10m, :pres, :dswrf, :dlwrf, :prate),
        area, years)

    # Normalise any 3-hourly layers (8/day) to 6-hourly (4/day).
    # See note in file header; SubDailyResolution{N} will handle this properly.
    n6h = sum(Dates.daysinyear, years) * 4
    names = keys(stack)
    layers = map(names) do name
        layer = getproperty(stack, name)
        n = size(layer, Ti)
        n == n6h && return layer
        n % n6h == 0 || error(
            "NCEP: $name has Ti=$n, not a multiple of expected 6h count $n6h")
        return _aggregate_ti_to_daily(layer, n ÷ n6h)
    end

    return Rasters.replace_missing(RasterStack(NamedTuple{names}(layers)), NaN)
end
