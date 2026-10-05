# Interactive

MicroclimateMapper.jl comes with a small interactive app that runs on your own computer, in your browser. Choose a
place on a map, swap the data set behind it, change a few settings, compare data sets, and copy a Julia script that
reproduces what you did, to build from.

```julia
using MicroclimateMapper, Bonito, WGLMakie
ENV["RASTERDATASOURCES_PATH"] = joinpath(homedir(), "spatial_data")   # where data sets are downloaded
server = MicroclimateMapper.app()        # opens http://127.0.0.1:8080 in your browser
# ...
close(server)
```

Start Julia with `--threads=auto` so that runs use every core. The first run of a session compiles the model and
takes a minute or two; later runs take seconds for a climatology, longer for a year of daily or hourly data.

## What it offers

| Control | What it changes |
|:--|:--|
| map and place search | the point, by clicking or by name ([`geocode`](@ref)); the data sets covering it are listed |
| data set | `weather_source`: CRU CL 2.0, WorldClim, TerraClimate, gridMET, SILO, NCEP, ERA5, ERA5-Land (`ERA5ECMWFLand`) |
| also run, to compare | further data sets at the same point, overlaid |
| year | the year, for data sets that are not climatologies |
| shade | `data = (; shade)` |
| snow | `SnowModel()` or `NoSnow()`; off by default |
| dynamic soil moisture | soil moisture solved from rainfall (`DynamicSoilMoisture()`), or prescribed from `CPCSoil`; off by default |
| terrain | off by default: flat, at the elevation of CRUCL2's 10′ cell, with slope, aspect and horizon from sliders; on: elevation, slope, aspect and horizon from SRTM (`dem_source`, `compute_terrain`) |
| albedo, soil | `surface_albedo_source`; the example loam, or a profile from SoilGrids or SLGA ([`build_soil_profile`](@ref)) |
| data folder, SILO email, Copernicus key | `ENV["RASTERDATASOURCES_PATH"]`, `ENV["SILO_EMAIL"]`, and point-by-point ERA5 from ECMWF's store, see [Access](manual/forcing_data.md#Access); kept in the session, never written to the script |

The model always has the same depths (surface to 2 m) and heights (1 cm to 2 m, the reference height). After a
run, the plots can be changed without running again:

- **what to plot**: soil temperature; air temperature, wind speed and humidity; radiation (solar, on a horizontal
  surface and on the slope, and longwave down and up); soil moisture, humidity and water potential (with dynamic
  soil moisture); snow (with snow); and dew and frost;
- **which depths and heights**, by checkbox;
- **the period**: the whole year, a month or a day, with sunrise and sunset marked on a single day;
- **the terrain around the point**: maps of elevation, slope, aspect and sky view factor from SRTM, and the
  horizon in every direction with the sun's path over it at the solstices and equinox, or on the day shown.

Results can be saved as CSV files, one per data set, in a `microclimate_outputs` folder in the data folder.

## The script

Every run writes a script with the same settings: the point, the data sets, the model and the problem. Copy it into
a file and run it to get the same `RasterStack`s, then go further than the app does: other depths and heights, a
grid instead of a point, other output layers, other years. See [Get started](get_started.md) and the
[tutorials](tutorials/points.md).

## Data and credentials

Data sets are downloaded the first time they are used, which for some, TerraClimate and gridMET above all, means
files of hundreds of megabytes; see [Forcing data sets](manual/forcing_data.md). `SILO` needs an email address, and `ERA5`
needs `using ZarrDatasets`; with a key for the Climate Data Store it is read point by point from ECMWF's store,
which is much faster. A data set that
fails shows its error in the app, and the others still run.

The app is the third of a set: SolarRadiation.jl and Microclimate.jl are to have interactive pages of their own,
computed in the browser, building up to this one.
