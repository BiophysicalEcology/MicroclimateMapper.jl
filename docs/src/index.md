# MicroclimateMapper.jl

Spatial driver for [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl).

MicroclimateMapper.jl runs the Microclimate microclimate model over rasters and
point sets, sourcing terrain (`SRTM`), weather (`TerraClimate`, `NCEP`, `ERA5`,
`WorldClim`, `CHELSA`, `AWAP`, `GRIDMET`), and land-cover data automatically.

## Problem types

- [`MicroRasterProblem`](@ref) — solve over a rectangular area at the
  resolution of a template raster.
- [`MicroVectorProblem`](@ref) — solve at a list of points.

Both share a [`MicroMapModel`](@ref) that wraps an inner Microclimate
`MicroModel` and the data sources to draw on.

## Quick start

```julia
using MicroclimateMapper
using Microclimate: example_microclimate_problem, example_soil_profile
using Dates

micro_model = example_microclimate_problem().model

model = MicroMapModel(micro_model;
    dem_source     = SRTM,
    weather_source = TerraClimate{Historical},
)

problem = MicroVectorProblem(model;
    points       = [(3.581, 44.122)],          # (lon, lat)
    dates        = Date(2000, 6, 1):Day(1):Date(2000, 6, 7),
    soil_profile = example_soil_profile(inner.depths),
)

output = solve(problem)
```

See the [Raster](examples/raster.md) and [Vector](examples/vector.md)
examples for larger runs, and the [API](api.md) for every exported type and
function.
