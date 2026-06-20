# Vector example

Point microclimate at three Cévennes sites that share the same NCEP grid
cell (~2.5°). Terrain lapse-correction, aspect, and elevation-driven
snowpack explain the differences between sites.

This example mirrors `vector_demo.jl` in
[MicroclimateTalk](https://github.com/BiophysicalEcology/MicroclimateTalk),
shortened to a winter window so it stays fast for documentation builds.

<!-- @setup vector
include(joinpath(@__DIR__, "..", "plotting.jl"))
-->

```julia
using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_profile,
    example_soil_properties_model, example_soil_hydraulic_model
using RasterDataSources
using Dates
using Unitful

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 2.0]u"m"

# Valley floor → south-facing slope → summit. Same NCEP cell, very
# different microclimates from elevation + aspect.
points = [
    (3.520, 44.143),   # Camprieu (valley, ~900 m)
    (3.582, 44.115),   # South-facing slope (~1200 m)
    (3.581, 44.122),   # Mont Aigoual summit (~1567 m)
]

dates = Date(2010, 2, 1):Day(1):Date(2010, 2, 28)

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths,
        heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(),
        snow_model            = SnowModel(),
    ),
    dem_source              = SRTM,
    weather_source          = NCEP,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
)

problem = MicroVectorProblem(;
    model,
    points,
    dates,
    soil_profile = example_soil_profile(depths),
    init         = (;
        soil_moisture = fill(0.25, length(depths)),
        snow_depth    = 0.0u"cm",
    ),
)

output = solve(problem)
```

The result is a `RasterStack` whose `point` dimension is a `MergedLookup`
of the input `(lon, lat)` tuples. Slicing along `point`, `Ti`, and (for
soil/air variables) `depth` / `height` gives per-site time series.

## Snow depth

One line per site, noon samples across the run window.

```julia
plot_vector_snow_depth(output,
    ["Valley (900 m)", "South slope (1200 m)", "Summit (1567 m)"])
```
