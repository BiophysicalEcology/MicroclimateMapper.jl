# ERA5 bindings: hourly reanalysis, for every RasterDataSources provider and store.
#
# - `ERA5` / `ERA5{Google,ARCO}`: Google's public Zarr store, no account needed.
# - `ERA5{ECMWF,ARCO}` / `ERA5{ECMWF,ARCO,Land}`: ECMWF's Zarr stores, no queue.
# - `ERA5{ECMWF,OnDemand}` / `ERA5{ECMWF,OnDemand,Land}`: monthly NetCDF
#   downloads from the Climate Data Store, for the area of the run.
#
# ECMWF needs a CDS token in `ENV["RASTERDATASOURCES_ERA5_CDS_TOKEN"]`.
# ARCO stores need `using Zarr`.
#
# Layers are the CDS names RasterDataSources uses for every store.

const _ERA5Land = ERA5{<:RasterDataSources.ERA5Provider,<:RasterDataSources.ERA5Store,RasterDataSources.Land}

weather_calendar(::Type{<:ERA5}) = Daily()
native_timestep(::Type{<:ERA5}) = Hourly()
loader(::Type{<:ERA5}) = ARCOTimeSeries()
loader(::Type{<:ERA5{ECMWF,OnDemand}}) = MonthlyTimeSeries()

# Google's store is 0-360°, ECMWF's stores and CDS downloads are -180-180°.
longitude_convention(::Type{<:ERA5}) = Longitude360()
longitude_convention(::Type{<:ERA5{ECMWF}}) = Longitude180()

# One cell around the points is enough: 0.25° for ERA5, 0.1° for ERA5-Land.
points_load_buffer(::Type{<:ERA5}) = 0.5
points_load_buffer(::Type{<:_ERA5Land}) = 0.2

# ECMWF files and stores name variables by their short names, e.g. `t2m`.
file_variable(T::Type{<:ERA5}, name::Symbol) = Symbol(layername(T, name))

_era5_variables() = (
    Variable(Reference(Temperature()), :temperature_2m, u"K"),
    Variable(EastwardWindSpeed(), :u_component_of_wind_10m, u"m/s"),
    Variable(NorthwardWindSpeed(), :v_component_of_wind_10m, u"m/s"),
    Variable(DewpointTemperature(), :dewpoint_temperature_2m, u"K"),
    Variable(Pressure(), :surface_pressure, u"Pa"),
    Variable(GlobalRadiation(), :surface_solar_radiation_downwards, u"J/m^2/hr"),
    Variable(LongwaveRadiation(), :surface_thermal_radiation_downwards, u"J/m^2/hr"),
    # Metres of water per hour; × 1000 kg/m³ gives kg/m².
    Variable(Rainfall(), :total_precipitation, u"kg/m^2", raw -> raw * 1000.0),
)

# Total cloud cover is a 0-1 fraction.
variables(::Type{<:ERA5}) = (_era5_variables()..., Variable(CloudCover(), :total_cloud_cover, 1))
# ERA5-Land has no cloud cover: it comes from ERA5 in the same store.
variables(::Type{<:_ERA5Land}) = _era5_variables()
# ERA5-Land downloads accumulate fluxes through each day (ECMWF's ARCO store
# is already hourly).
const _ERA5_LAND_ACCUMULATED = (:surface_solar_radiation_downwards,
    :surface_thermal_radiation_downwards, :total_precipitation)
_from_accumulated(::Type{<:ERA5{ECMWF,OnDemand,RasterDataSources.Land}}, name::Symbol, layer) =
    name in _ERA5_LAND_ACCUMULATED ? _hourly_from_daily_accumulation(layer) : layer

fallback_source(::Type{<:ERA5{ECMWF,S,RasterDataSources.Land}}) where {S<:RasterDataSources.ERA5Store} = ERA5{ECMWF,S}
