# ERA5 bindings — hourly global reanalysis from the ARCO-ERA5 cloud-optimised
# Zarr store. 0.25° resolution, 1940-present (~3 month lag), 8760 timesteps
# per year (leap-day Feb 29 dropped).
#
# In hourly mode we bypass the min/max envelope (`environment_minmax = nothing`)
# and let the solver dispatch on the `::Nothing` path, which reads the five
# required fields straight from `environment_hourly`. The canonical buffers
# `reference_temperature`, `reference_humidity`, `wind_speed`, `cloud_cover`,
# `pressure`, `global_radiation`, `longwave_radiation`, `rainfall`,
# `zenith_angle` are shared storage with `environment_hourly`'s fields, so
# writes via the canonical name surface in the struct.
#
# The hourly chain (`_DERIVATIONS_HOURLY`) handles four conversions:
#   * `wind_speed`              ← √(u_wind² + v_wind²)   (from `:u10`/`:v10`)
#   * `actual_vapour_pressure`  ← e_s(dewpoint_temperature) (from `:d2m`)
#   * `reference_humidity`      ← actual_VP / e_s(reference_temperature)
#   * `rainfall_daily_from_hourly` and `deep_soil_temperature_from_hourly`
#     aggregate hourly canonicals into the daily `environment_daily` fields.

temporal_resolution(::Type{<:ERA5}) = HourlyResolution()
weather_loader(::Type{<:ERA5})      = HourlyZarrStore()
weather_derivations(::Type{<:ERA5}) = _DERIVATIONS_HOURLY

function weather_variables(::Type{<:ERA5})
    (
        WeatherVariable(:reference_temperature, :t2m,  u"K"),
        WeatherVariable(:u_wind,                :u10,  u"m/s"),
        WeatherVariable(:v_wind,                :v10,  u"m/s"),
        WeatherVariable(:dewpoint_temperature,  :d2m,  u"K"),
        WeatherVariable(:pressure,              :sp,   u"Pa"),
        # Total cloud cover already 0-1 fraction; no unit, no transform.
        WeatherVariable(:cloud_cover,           :tcc,  1),
        # ERA5 radiation is stored as J/(m²·hour) accumulations.
        # Unitful converts that to W/m² on assignment to the canonical buffer.
        WeatherVariable(:global_radiation,      :ssrd, u"J/m^2/hr"),
        WeatherVariable(:longwave_radiation,    :strd, u"J/m^2/hr"),
        # `:tp` is total precipitation in metres of water per hour; multiply
        # by 1000 to get kg/m² (assuming water density of 1000 kg/m³).
        WeatherVariable(:rainfall, :tp, u"kg/m^2", raw -> raw * 1000.0),
    )
end
