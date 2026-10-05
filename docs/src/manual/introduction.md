# Introduction

Microclimate modellers want small-scale predictions over large areas, often hourly. That is a lot of
computation, and it competes with the need to try new ideas quickly. The answer here is a set of *levers* that
switch between research modes without rewriting anything:

- **spatial resolution**: the grid a run is solved on;
- **temporal resolution**: monthly, daily or hourly forcing;
- **formulation accuracy**: which equation for each process;
- **formulation completeness**: which processes, such as snow or a canopy, are included.

[Microclimate.jl](https://biophysicalecology.github.io/Microclimate.jl/dev) solves the microclimate at a point
from weather, terrain, soil and vegetation inputs. MicroclimateMapper.jl supplies those inputs from spatial data
sets, for one point, many points or every cell of a grid.

## Model and problem

The packages of [BiophysicalEcology](https://github.com/BiophysicalEcology) use two words throughout:

- a **model** is a physical specification, without a place or time;
- a **problem** is a model with a place and time, which can be solved.

A [`MicroMapModel`](@ref) wraps the point model of Microclimate.jl, a `MicroModel`, with the data sets to drive
it: a DEM, a weather source, and optionally land cover, surface properties and soil moisture. It says nothing
about where or when. A [`MicroVectorProblem`](@ref) pairs it with points and dates; a
[`MicroRasterProblem`](@ref) with an area, a grid and dates. `solve` returns a `RasterStack`.

```julia
model = MicroMapModel(; micro_model, dem_source = SRTM, weather_source = CRUCL2)
points = MicroVectorProblem(; model, points = [(-89.40, 43.07)], dates, soil_profile)
grid = MicroRasterProblem(; model, area, template = SRTM, dates, soil_profile)
solve(points), solve(grid)
```

The design follows from a few principles:

- **Declarative**: a script says what to solve, not how.
- **Modular**: data sets and model components are swapped for alternatives with one keyword, see
  [One keyword per data set](data_sources.md).
- **Unit-safe**: every quantity carries its units, which are checked.
- **Readable**: formulations are short and look like the equations they come from.
- **Fast**: Julia can approach the performance of Fortran and C, see [Performance](performance.md).

## How it combines with the Microclimate.jl package

| Step | Done here |
|:--|:--|
| fetching data | data sets as types of [RasterDataSources.jl](https://ecojulia.github.io/RasterDataSources.jl/stable), downloaded once and cached, see [Built on Rasters and DimensionalData](spatial_stack.md) |
| forcing | each data set's variables converted to one hourly forcing, whatever its native time step, see [Forcing data sets](forcing_data.md) |
| terrain | slope, aspect and horizon angles from any DEM, with [Geomorphometry.jl](https://deltares.github.io/Geomorphometry.jl/stable/), and weather corrected for elevation, see [Terrain and atmosphere](terrain.md) |
| surfaces and soils | albedo and roughness from land cover, soil profiles from texture maps, see [Surfaces and soils](surfaces_soils.md) |
| running | every point or cell in parallel, with the results gathered into one `RasterStack`, see [Outputs](outputs.md) |

Microclimate.jl treats each point as part of an infinite uniform plane. Flows between neighbouring cells, such
as runoff and cold-air drainage, need the grid and so belong here. Runoff routing is in development, see
[Lateral flows](lateral_flows.md).

## The ecosystem

| Package | Role |
|:--|:--|
| [Microclimate.jl](https://biophysicalecology.github.io/Microclimate.jl/dev) | the microclimate at a point |
| [SolarRadiation.jl](https://biophysicalecology.github.io/SolarRadiation.jl/dev) | clear-sky solar radiation over terrain |
| [FluidProperties.jl](https://biophysicalecology.github.io/FluidProperties.jl/stable/) | properties of air and water |
| [HeatExchange.jl](https://biophysicalecology.github.io/HeatExchange.jl/dev) | heat and water budgets of organisms in these microclimates |
| [BiophysicalBehaviour.jl](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev) | organisms choosing among these microclimates, see [Activity and available environments](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/environments#From-gridded-climate-data) |
| [RasterDataSources.jl](https://ecojulia.github.io/RasterDataSources.jl/stable), [Rasters.jl](https://rafaqz.github.io/Rasters.jl/stable), [DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl/stable/), [Geomorphometry.jl](https://deltares.github.io/Geomorphometry.jl/stable/), [PointDataSources.jl](https://github.com/BiophysicalEcology/PointDataSources.jl) | fetching, holding and analysing the spatial data |

MicroclimateMapper.jl is the Julia counterpart of the `micro_global`, `micro_usa`, `micro_ncep` and related
functions of [NicheMapR](https://github.com/mrke/NicheMapR) (Kearney and Porter 2017), see
[For NicheMapR users](nichemapr.md).