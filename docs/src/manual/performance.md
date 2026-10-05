# Performance

The levers of the [Introduction](introduction.md) are also levers of cost. A run solves the point model once for
every point or cell, for every day, and holds every output layer for every hour in memory. Its cost is roughly

```math
\text{cells} \times \text{days} \times \text{cost per cell-day}
```

where the days are 12 for a monthly data set and 365 a year for a daily or hourly one.

## What it costs

Times measured on a laptop in October 2026 (Windows, 16 threads):

| Run | Cells | Days | Solve | Per cell-day |
|:--|:--|:--|:--|:--|
| two points, CRU CL 2.0, with snow | 2 | 12 | 0.6 s | 26 ms (one thread per point) |
| Palm Springs–San Jacinto, ~900 m, CRU CL 2.0, with snow | 2268 | 12 | 284 s | 10 ms |
| the same, without snow | 2268 | 12 | 185 s | 7 ms |
| California at 10′, CRU CL 2.0, with snow | 3591 | 12 | 456 s | 11 ms |

So on this machine, a grid of a few thousand cells for a monthly climatology takes minutes, and the same grid for a
year of daily weather takes hours.

What raises the cost of each cell-day:

- **snow**, by about half again;
- **dynamic soil moisture**, which solves infiltration and redistribution every hour;
- **a multilayer canopy**, which iterates the energy balance of every leaf layer;
- **more depths and heights**.

## Memory

The output of a grid is held in memory: cells × hours × (depths or heights) for each layer. A 100 × 100 grid of a
monthly climatology, 288 hours with the default layers, has run out of memory on a laptop. Write only the layers
needed, through `output_layers`, see [Outputs](outputs.md). Split a large area into tiles and write each to disk.

## Loading

`init` loads and prepares the data, and can take as long as the solve for a small run:

- **The first download** of a data set is the slowest step of all. Later runs read the cache.
- **Terrain**: horizon angles are computed on the DEM at its own resolution, which is costly for a fine DEM over a
  large area. Reuse it with `data = (; terrain = terrain(cache))`, or skip it with `compute_terrain = false` on
  coarse grids.
- **Points far apart** share one DEM over the box that contains them, see [Points](../tutorials/points.md#Points-far-apart).

## Running again

`solve(problem)` is `solve!(init(problem))`. Keeping the cache separates loading from solving, and `solve!` can be
called again on the same cache:

```julia
cache = init(problem)
output = solve!(cache)
```

Each thread has its own point-model cache, reset for each cell rather than rebuilt. Start Julia with
`--threads=auto`, or set `JULIA_NUM_THREADS`, to use every core.

## Fetching ahead

[`prefetch_weather!`](@ref)`(source, points, dates)` downloads a data set for the points' area and dates in
batches of a week, skipping batches already cached, so an interrupted download resumes. It suits hourly data sets
read from servers, such as `ERA5` and `BARRA`. Read directly from Google's ARCO-ERA5 store, which holds the whole
globe one hour at a time, a year of hourly `ERA5` at one point took more than 25 minutes in the run behind
[Swapping data sets](../tutorials/swap_data.md). With a Copernicus API key, `ERA5ECMWF` reads ECMWF's own stores instead,
and for points ECMWF's copy chunked by place, through the point queries below, is much faster still, see
[Forcing data sets](forcing_data.md#Access).

## Point queries

For a few points, downloading grids is wasteful. The `PointQuery` loader fetches each variable's time series at
each point through [PointDataSources.jl](https://github.com/BiophysicalEcology/PointDataSources.jl) instead. It
is chosen per data set:

```julia
MicroclimateMapper.loader(::Type{<:SILO}) = MicroclimateMapper.PointQuery()
MicroclimateMapper.loader(::Type{<:ERA5}) = MicroclimateMapper.PointQuery()   # with ENV["CDS_API_KEY"] set
```

The [interactive app](../interactive.md) reads ERA5 this way when given a key.

## Inspecting before solving

With Makie loaded, `plot(cache)` shows the inputs that `init` has prepared, weather, terrain and surface, with a
slider through time, before any solving is done.
