# Built on Rasters and DimensionalData

MicroclimateMapper.jl does little spatial work of its own. It stands on four packages of the
[JuliaGeo](https://juliageo.org/) and [EcoJulia](https://github.com/EcoJulia) ecosystems:

| Package | Does here |
|:--|:--|
| [DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl/stable/) | arrays whose dimensions have names and coordinates |
| [Rasters.jl](https://rafaqz.github.io/Rasters.jl/stable) | rasters and stacks of them: reading, cropping, resampling, masking, writing |
| [RasterDataSources.jl](https://ecojulia.github.io/RasterDataSources.jl/stable) | data sets as types, downloaded once and cached |
| [Geomorphometry.jl](https://deltares.github.io/Geomorphometry.jl/stable/) | slope, aspect and horizon angles from elevation |

and on [PointDataSources.jl](https://github.com/BiophysicalEcology/PointDataSources.jl) for time series at points.
This page shows what each contributes, on the area around Palm Springs of [Get started](../get_started.md).

```@setup stack
using Main.FigureHelpers
using CairoMakie
```

## Data sets as types

A data set is asked for by its type, not by a URL or a file format. The file is downloaded the first time, into
`ENV["RASTERDATASOURCES_PATH"]`, and read from there afterwards:

```@example stack
using Rasters, RasterDataSources
using Rasters.Extents: Extent
area = Extent(X = (-116.80, -116.45), Y = (33.70, 33.95))
srtm = Raster(SRTM; extent = area, lazy = true)
```

`lazy = true` reads nothing until the raster is indexed, so a tile or a global file costs no memory until a part
of it is used. Other data sets work the same way: `Raster(CHELSA{Climate}, :tas; month = 7)`,
`RasterStack(TerraClimate{Historical}, (:tmax, :tmin); date)`, `RasterStack(NCEP{SurfaceFlux, 1}; date)`. In R,
each would be a different function or package, with its own arguments and file handling.

## Named dimensions

A raster is an array whose dimensions have names and coordinates. It is indexed by what is meant, not by
position:

```@example stack
srtm[X(Near(-116.545)), Y(Near(33.830))]    # Palm Springs
```

and a region is selected by coordinate range:

```@example stack
peak = srtm[X(-116.75 .. -116.60), Y(33.75 .. 33.88)]
size(peak)
```

In R's terra the same is `extract(srtm, cbind(-116.545, 33.830))` and `crop(srtm, ext(...))`; with an array,
`elev[150, 200]`, and the row and column must be worked out by hand. Here a mistake of dimension is an error,
much as [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) makes a mistake of units an error.

The outputs of this package are the same kind of object, with `point` or `X` and `Y`, `Ti`, and `depth` or
`height`:

```julia
output.soil_temperature[point = 1, depth = Near(0.1)]   # a time series
output.air_temperature[Ti = Near(DateTime(2000, 7, 1, 14))]  # a map, or a value per point
maximum(output.soil_temperature; dims = Ti)                  # the hottest hour at every depth and place
```

Any tool built on DimensionalData.jl can read them, and [Makie](https://docs.makie.org/) plots them directly.

## Terrain

[Geomorphometry.jl](https://deltares.github.io/Geomorphometry.jl/stable/) computes slope and aspect from the eight
neighbours of each cell (Horn 1981), and the angle of the horizon in each direction. It accepts rasters, and
reads their coordinate reference system to get the cell size in metres:

```@example stack
using Geomorphometry
using Geomorphometry: Horn
dem = Float64.(read(peak))                      # Geomorphometry.jl works on real numbers
slope_angle = Geomorphometry.slope(dem; method = Horn())
slope_angle = slope_angle[X(2:size(dem, X) - 1), Y(2:size(dem, Y) - 1)]   # the edge cells lack neighbours
fig = Figure(size = (800, 330)) # hide
map_panel!(fig[1, 1], dem; title = "Elevation (m)", colormap = :terrain, label = "m") # hide
map_panel!(fig[1, 2], slope_angle; title = "Slope (°)", colormap = :viridis, label = "°") # hide
fig # hide
```

This is done for every run, from the run's `dem_source`, see [Terrain and atmosphere](terrain.md). SolarRadiation.jl's
[Terrain: Saba](https://biophysicalecology.github.io/SolarRadiation.jl/dev/tutorials/saba) does the same by hand,
step by step.

## Grids and resampling

A [`MicroRasterProblem`](@ref) solves on a `template` grid, and everything else is brought onto it: the weather
by resampling from its own grid, land cover and surface properties likewise, and the terrain by aggregating from
the finer DEM. A coarser grid is one call:

```@example stack
using Statistics
coarse = aggregate(mean, dem, 10)
size(dem), size(coarse)
```

which is the *spatial resolution* lever of the [Introduction](introduction.md): the same run on `coarse` costs a
hundredth as much.

An `area` can also be any polygon, from a shapefile or GeoJSON. Cells outside it are skipped, and left `missing`
in the output.

## Writing

Outputs carry units, which NetCDF cannot store. [`strip_to_canonical`](@ref) converts each layer to a standard
unit and strips it, so that Rasters.jl can write the stack:

```julia
write("palm_springs.nc", strip_to_canonical(output))
```

## Points

For a few points, [PointDataSources.jl](https://github.com/BiophysicalEcology/PointDataSources.jl) fetches time
series from a data set's point service, `getpoint(SILO, :max_temp; lon, lat, date)`, instead of downloading grids.
The `PointQuery` loader makes a data set use it, see [Performance](performance.md#Point-queries).
