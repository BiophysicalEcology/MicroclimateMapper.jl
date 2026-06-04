# Raster example

Gridded microclimate over a small high-resolution box centred on
Mont Aigoual (Cévennes, southern France) at the native SRTM resolution,
forced by TerraClimate.

This example mirrors `raster_demo.jl` in
[MicroclimateTalk](https://github.com/BiophysicalEcology/MicroclimateTalk),
shortened to a one-week window so it stays fast for documentation builds.

```@setup raster
include(joinpath(@__DIR__, "..", "plotting.jl"))
```

```@example raster
using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_profile,
    example_soil_properties_model, example_soil_hydraulic_model
using RasterDataSources
using Rasters
using Rasters.Extents: Extent
using Dates
using Unitful

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 2.0]u"m"

# Mont Aigoual: 1567 m peak, dramatic N/S aspect contrast.
# ~0.054° square centred on the summit.
area  = Extent(X = (3.554, 3.608), Y = (44.095, 44.149))
dates = Date(2000, 7, 1):Day(1):Date(2000, 7, 7)

micro_model = MicroModel(;
    depths,
    heights,
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model  = example_soil_hydraulic_model(),
    snow_model            = NoSnow(),
)

map_model = MicroMapModel(;
    micro_model,
    dem_source              = SRTM,
    weather_source          = TerraClimate{Historical},
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
)

problem = MicroRasterProblem(;
    model        = map_model,
    area,
    dates,
    template     = SRTM,
    soil_profile = example_soil_profile(depths),
)

output = solve(problem)
nothing # hide
```

The result is a `RasterStack` keyed by output layer, with `X`, `Y`, `Ti`,
and (where applicable) `depth` / `height` dimensions. For longer runs the
talk's demo extends `dates` over a full year and writes to NetCDF.

## Soil temperature at 5 cm

Six panels show the same patch at midnight, dawn, mid-morning, midday,
mid-afternoon, and dusk on the first day of the run, draped on the SRTM
DEM. Aspect and elevation drive the spread across cells.

```@example raster
soil_temperature_5cm = ustrip.(u"°C",
    view(output.soil_temperature; depth = 3))

dem = Rasters.resample(
    read(crop(Raster(SRTM; extent = area, lazy = true, missingval = Int16(0));
              to = area, touches = true));
    to = first(output),
)

plot_raster_temperature_panels(soil_temperature_5cm, dem)
```
