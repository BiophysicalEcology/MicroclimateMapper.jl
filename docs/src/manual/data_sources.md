# One keyword per data set

Changing the weather that drives a run is one keyword:

```julia
model = MicroMapModel(; micro_model, dem_source = SRTM, weather_source = CRUCL2)
model = MicroMapModel(; micro_model, dem_source = SRTM, weather_source = NCEP{SurfaceFlux, 1})
model = MicroMapModel(; micro_model, dem_source = SRTM, weather_source = GRIDMET)
```

Nothing else needs to change to produce either a monthly climatology, a 6-hourly reanalysis or a daily gridded
record, with different variables, units, grids and file layouts. This page explains how, and how to add a data
set.

```@setup sources
using Main.FigureHelpers
```

## Multiple dispatch

Julia chooses which method of a function to run from the types of *all* its arguments. Every data set is a type,
so one function can have a method per data set. How the files of a data set are laid out, for example:

```@example sources
using MicroclimateMapper, RasterDataSources
const M = MicroclimateMapper
M.loader(CRUCL2), M.loader(GRIDMET), M.loader(ERA5)
```

CRUCL2 is one file of monthly bands, GRIDMET a file per year and ERA5 one contiguous cloud store. The code that
reads the weather calls `loader(source)` and goes the right way for each. A method you define
takes precedence over the package's, so how a data set is read can be changed without touching the package. By
default ERA5 is read from Google's ARCO copy, a cloud store chunked by time, which is slow for a single point.
With a Copernicus API key, one line reads it point by point from ECMWF's own store, chunked by place (this is what
the [interactive app](../interactive.md) does):

```julia
using ZarrDatasets
ENV["CDS_API_KEY"] = "<your key>"
MicroclimateMapper.loader(::Type{<:ERA5}) = MicroclimateMapper.PointQuery()
```

Every model using `weather_source = ERA5` then reads it that way. The same pattern runs through the BiophysicalEcology
packages, for formulations as well as data: see FluidProperties.jl's
[Adding an equation](https://biophysicalecology.github.io/FluidProperties.jl/dev/manual/vapour_pressure#Adding-an-equation)
for a vapour pressure equation of your own.

## Data sets are types

[RasterDataSources.jl](https://ecojulia.github.io/RasterDataSources.jl/stable) gives every data set a type:
`CRUCL2`, `TerraClimate`, `GRIDMET`, `NCEP{SurfaceFlux, 1}`, `ERA5` and so on. `Raster(GRIDMET, :tmmx; date)`
downloads the file if it is not already cached and returns the layer as a raster. MicroclimateMapper.jl adds a
few functions of the data-set type, *traits*, that say what the model needs to know about it:

| Trait | Says | Values |
|:--|:--|:--|
| `weather_calendar` | which days the data describe | `Monthly()` (one representative day per month) or `Daily()` |
| `native_timestep` | the samples within a day | `MinMax()` (daily extremes), `SixHourly()` or `Hourly()`; these two are `SubDaily{4}()` and `SubDaily{24}()`, and any `SubDaily{N}()` (N samples a day) works |
| `loader` | how the files are laid out | one file per year, per month or per day, one file of bands, a contiguous cloud store or several (one per group of variables), point queries |
| `variables` | which layers to read, as what, in which units | below |
| `fallback_source`, `fallback_layers` | where to get what the data set lacks | e.g. wind for `AWAP` and `SILO` from `CRUCL2` |

For the data sets available now:

```@example sources
timestep(t) = t isa M.MinMax ? "`MinMax()`, daily extremes" : t isa M.Hourly ? "`Hourly()`" :
    t isa M.SixHourly ? "`SixHourly()`" : "`SubDaily{$(M.samples_per_day(t))}()`"
fallback(S) = (f = M.fallback_source(S)) === nothing ? "—" : string(f)
sources = (CRUCL2, WorldClim{Climate}, CHELSA{Climate}, TerraClimate, GRIDMET, NCEP{SurfaceFlux, 1}, ERA5, ERA5ECMWF,
    ERA5ECMWFLand, AWAP, SILO, BARRA)
markdown_table(["Data set", "Calendar", "Within a day", "Loader", "Fallback"],
    [("`$S`", nameof(typeof(M.weather_calendar(S))), timestep(M.native_timestep(S)),
      nameof(typeof(M.loader(S))), fallback(S)) for S in sources])
```

## Canonical variables

Each data set declares its layers as physical quantities, with the name and unit of the layer in its files.
TerraClimate's:

```@example sources
markdown_table(["Quantity", "Canonical name", "Layer in the files", "Unit"],
    [(nameof(typeof(M.quantity(v))), "`$(M.canonical_name(v))`", "`$(M.native_field(v))`", v.unit)
     for v in M.variables(TerraClimate)])
```

A quantity is a type such as `Temperature(Maximum())`, `RelativeHumidity(Minimum())`,
`ActualVapourPressure(ClockTime(9))` or `Reference(Temperature())`, the last meaning measured at the reference
height. Its canonical name is the same whichever data set it came from. Units are converted on reading, and a
`transform` handles anything else: TerraClimate's soil water in millimetres becomes a volumetric fraction, NCEP's
precipitation rate becomes a 6-hourly total.

The model needs a fixed set of forcing variables. What a data set lacks is derived from what it has, by methods
of one function, `derive!`, chosen by dispatch on the quantity:

- extremes from a mean and a range (`CRUCL2` gives mean temperature and diurnal range);
- vapour pressure from relative humidity, vapour pressure deficit, dew point, specific humidity, or readings at
  09:00 and 15:00 (`AWAP`);
- wind speed from its eastward and northward components (`NCEP`, `ERA5`);
- daily extremes from sub-daily samples, and hourly values from daily extremes through the diel curves of
  Microclimate.jl, see its [Diel curves and time handling](https://biophysicalecology.github.io/Microclimate.jl/dev/manual/diel_time_handling).

So a data set declares only what it has, and the derivations are shared by all.

## Overriding

Any canonical variable can be supplied directly through `data` on a problem, as a `Raster` in canonical units,
in place of what the data set gives. It is resampled to the run's grid or points:

```julia
problem = MicroRasterProblem(; model, area, template, dates, soil_profile,
    data = (; cloud_cover = my_cloud_raster, shade = 0.3))
```

`data` also takes a whole weather stack (`weather`), a DEM (`dem`), terrain from an earlier run (`terrain`),
land cover, albedo and roughness, and the fraction of `shade`. Data supplied this way take precedence over the
data sets of the model.

## Adding a data set

A data set already in RasterDataSources.jl needs a calendar, a loader and its variables. TerraClimate's binding
is the whole of `src/climate/terraclimate.jl`:

```julia
loader(::Type{<:TerraClimate}) = YearlyTimeSeries()   # one file per year; calendar is Monthly() by default

function variables(::Type{<:TerraClimate})
    (
        Variable(Temperature(Maximum()), :tmax, u"°C"),
        Variable(Temperature(Minimum()), :tmin, u"°C"),
        Variable(WindSpeed(), :ws, u"m/s"),
        Variable(VapourPressureDeficit(), :vpd, u"kPa"),
        Variable(GlobalRadiation(), :srad, u"W/m^2"),
        Variable(Rainfall(), :ppt, u"kg/m^2"),
        Variable(SoilMoisture(), :soil, 1, _terraclimate_soil_to_volumetric),
    )
end
```

A data set not yet in RasterDataSources.jl is added there first, see its
[documentation](https://ecojulia.github.io/RasterDataSources.jl/stable). A quantity no data set has yet, and no
`derive!` method can make, needs a new `derive!` method.
