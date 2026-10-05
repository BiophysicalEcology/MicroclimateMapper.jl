# Forcing data sets

The weather that drives a run comes from one keyword, `weather_source`, set to a data-set type of
[RasterDataSources.jl](https://ecojulia.github.io/RasterDataSources.jl/stable). Each data set reaches the model
as the same hourly forcing, see [One keyword per data set](data_sources.md). They differ in coverage,
resolution, period, time step and the variables they provide:

| `weather_source` | Coverage | Resolution | Period | Time step | Variables | Filled from elsewhere | NicheMapR |
|:--|:--|:--|:--|:--|:--|:--|:--|
| `CRUCL2` | global land | 10′ (~18 km) | 1961–1990 climatology | monthly means | mean temperature and diurnal range, humidity, wind, rainfall, sunshine, elevation | — | `micro_global` |
| `WorldClim{Climate}` | global land | 30″–10′ | 1970–2000 climatology | monthly | maximum and minimum temperature, wind, radiation, rainfall, vapour pressure | — | — |
| `CHELSA{Climate}` | global land | 30″ (~1 km) | 1981–2010 climatology | monthly | maximum and minimum temperature, wind, radiation, rainfall, humidity, cloud | — | — |
| `TerraClimate{Historical}` | global land | 1/24° (~4 km) | 1958 onwards | monthly | maximum and minimum temperature, wind, vapour pressure deficit, radiation, rainfall, soil moisture | — | `micro_terra` |
| `GRIDMET` | contiguous USA | 1/24° (~4 km) | 1979 onwards | daily | maximum and minimum temperature and humidity, wind, radiation, rainfall | — | `micro_usa` |
| `NCEP` (surface fluxes) | global | ~1.9° | 1948 onwards | 6-hourly | temperature, wind components, specific humidity, pressure, shortwave and longwave radiation, rainfall | — | `micro_ncep` |
| `ERA5` | global | 0.25° (~30 km) | 1940 onwards | hourly | temperature, wind components, dew point, pressure, cloud, shortwave and longwave radiation, rainfall | — | `micro_era5` |
| `ERA5ECMWF` | global | 0.25° (~30 km) | 1940 onwards | hourly | as `ERA5`, from ECMWF's own store (needs a Copernicus key) | — | `micro_era5` |
| `ERA5ECMWFLand` | global land | 0.1° (~9 km) | 1950 onwards | hourly | as `ERA5` but cloud | cloud from `ERA5ECMWF` | — |
| `AWAP` | Australia | 0.05° (~5 km) | 1900 onwards | daily | maximum and minimum temperature, rainfall, radiation, vapour pressure at 09:00 and 15:00 | wind from `CRUCL2` | `micro_aust` |
| `SILO` | Australia | 0.05° (~5 km) | 1889 onwards | daily | maximum and minimum temperature, rainfall, humidity at the times of the extremes, radiation | wind from `CRUCL2` | `micro_silo` |
| `BARRA` | Australia | ~12 km (R2), ~4 km (C2) | 1979 onwards | hourly | temperature, wind, humidity, sea-level pressure, shortwave and longwave radiation, rainfall, elevation | — | `micro_barra` |

!!! warning "Status, October 2026"
    The bindings are being brought in line with new versions of RasterDataSources.jl and Rasters.jl. In the
    versions used to build this documentation, `NCEP` and `BARRA` fail on loading (a layer name and a static
    elevation layer), and a year of hourly `ERA5` at a point takes tens of minutes to read from Google's store
    (use [`prefetch_weather!`](@ref), or the faster point route under [Access](#Access)). `CRUCL2`, `WorldClim`, `TerraClimate` and `SILO` run in the examples here.

Future climates come in two forms. `TerraClimate{Plus2C}` and `TerraClimate{Plus4C}` are TerraClimate's monthly
records with the climate shifted to 2 °C and 4 °C of global warming, with all its variables, and run as
TerraClimate does. `WorldClim{Future{Climate, …}}` and `CHELSA{Future{Climate, …}}` are climatologies for
projected periods under particular climate models and scenarios; they provide only temperature and rainfall, and
take their other variables from the matching baseline climatology.

## Choosing one

- **For a climatology anywhere**, `CRUCL2` is one small download and runs fastest: one representative day per
  month. It is the forcing of most examples in this documentation. `WorldClim` and `CHELSA` are finer.
- **For future climates**, `TerraClimate{Plus2C}` or `TerraClimate{Plus4C}` for warming levels, or the future
  climatologies of `WorldClim` and `CHELSA` for particular climate models and scenarios.
- **For particular years**, `TerraClimate` is global and monthly. `GRIDMET` (USA), `AWAP` and `SILO` (Australia)
  are daily and finer.
- **For weather hour by hour**, `ERA5` or `NCEP` globally, `BARRA` for Australia. Hourly forcing carries the
  timing of fronts, cloud and rain, which daily or monthly forcing can only approximate with a typical diel
  cycle. The hourly forcing includes the natural covariation among variables which can be important for phenomena 
  like dew and frost; such covariation is stereotyped for the daily and monthly forcing (e.g., lowest wind speed 
  is always at the same time of day).
- **For snow**, daily or hourly forcing. A monthly climatology has no individual storms, so snow cover comes and
  goes smoothly with the seasons.

Monthly and climatological sources solve 12 representative days, one per month, whatever the dates asked for.
Daily and sub-daily sources solve every day, carrying soil and snow state from one day to the next. The second
costs about thirty times as much for a year, see [Performance](performance.md).

## Access

Most sources are downloaded once and cached under `ENV["RASTERDATASOURCES_PATH"]`, see
[Built on Rasters and DimensionalData](spatial_stack.md). Some need more:

- **ERA5** can be reached three ways through RasterDataSources.jl, some still in open pull requests. All
  need the Zarr extension, `using ZarrDatasets`, except the last.
  - `weather_source = ERA5` reads Google's public ARCO-ERA5 store: no key and nothing downloaded beforehand,
    but the store is chunked one hour of the whole globe at a time, so a long record at one place is slow.
  - `weather_source = ERA5ECMWF` (or `ERA5ECMWFLand`, the 9 km land product) reads ECMWF's own analysis-ready
    stores, with no queue, given a Copernicus Climate Data Store API key in `ENV["CDS_API_KEY"]`. ECMWF calls
    this access a beta service. For points there is also a copy chunked by place, which is fastest for a long
    record at one point: the `PointQuery` loader reads it through PointDataSources.jl, with
    `MicroclimateMapper.loader(::Type{<:ERA5}) = MicroclimateMapper.PointQuery()`, see
    [One keyword per data set](data_sources.md#Multiple-dispatch).
  - `ERA5CDS` (and `ERA5CDSLand`) use the Climate Data Store's own API, which queues a request and downloads
    regional monthly files, as NicheMapR's `micro_era5` does through mcera5. Also keyed. Not yet bound to the
    model.
- **`SILO`** asks for an email address, read from `ENV["SILO_EMAIL"]`.
- **`BARRA`** is read from the THREDDS server of Australia's National Computational Infrastructure.
- **Point queries.** For a few points, the `PointQuery` loader fetches time series through
  [PointDataSources.jl](https://github.com/BiophysicalEcology/PointDataSources.jl) instead of grids, see
  [Performance](performance.md#Point-queries).

## Other sources

| Keyword | Data sets | Role |
|:--|:--|:--|
| `dem_source` | `SRTM` (~90 m), `CopernicusDEM` (30 m), `CRUCL2`, `BARRA` (their own elevation grids) | terrain and the elevation correction of the weather, see [Terrain and atmosphere](terrain.md) |
| `soil_moisture_source` | `CPCSoil` (monthly climatology, 0.5°) | prescribed soil moisture where the weather source has none, see [Surfaces and soils](surfaces_soils.md) |
| `init_source` | `ERA5` | initial soil temperature and moisture from one snapshot at the start date |
| `landcover_source` | `EarthEnv{LandCover}`, MODIS | albedo and roughness by land-cover class |
