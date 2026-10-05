# Get started

MicroclimateMapper.jl computes microclimates with [Microclimate.jl](https://biophysicalecology.github.io/Microclimate.jl/dev)
for points or grids, taking the weather, terrain, vegetation, soil and surface from spatial data sets.

```julia
using Pkg
Pkg.add(url = "https://github.com/BiophysicalEcology/MicroclimateMapper.jl")
```

Data sets are downloaded once and kept in a folder of your choice:

```julia
ENV["RASTERDATASOURCES_PATH"] = joinpath(homedir(), "spatial_data")
```

```@setup get_started
using Main.FigureHelpers
using CairoMakie
```

## A desert and a mountain

Palm Springs, California, lies at the foot of Mount San Jacinto, 13 km away and 3 km higher. The desert iguana of
[A lizard's day](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/tutorials/lizard) lives on the
desert floor.

The point model is a `MicroModel` of Microclimate.jl: soil depths, heights above ground, soil properties and
snow, see its [Configuring the model](https://biophysicalecology.github.io/Microclimate.jl/dev/manual/configuring_the_model).

```@example get_started
using MicroclimateMapper, Microclimate, RasterDataSources, Rasters, Dates, Unitful
using Microclimate: example_soil_profile, example_soil_properties_model, example_soil_hydraulic_model

depths = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"
micro_model = MicroModel(; depths, heights,
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
    snow_model = SnowModel(),
)
nothing # hide
```

A [`MicroMapModel`](@ref) adds the data sets. Here the weather is the CRU CL 2.0 climatology, soil moisture is
the CPC climatology, and terrain comes from SRTM:

```@example get_started
model = MicroMapModel(; micro_model,
    dem_source = SRTM,
    weather_source = CRUCL2,
    soil_moisture_source = CPCSoil,
    surface_albedo_source = 0.15,
    roughness_height_source = 0.004u"m",
)
nothing # hide
```

A [`MicroVectorProblem`](@ref) gives it places and dates. Points are `(longitude, latitude)`:

```@example get_started
sites = [(-116.545, 33.830), (-116.679, 33.814)]   # Palm Springs, San Jacinto Peak
problem = MicroVectorProblem(; model, points = sites,
    dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31),
    soil_profile = example_soil_profile(depths),
    init = (; soil_moisture = fill(0.2, length(depths)), snow_depth = 0.0u"cm"),
)
cache = init(problem)     # loads the data
output = solve!(cache)    # solves every point
```

`solve(problem)` does both. The output is a `RasterStack`. Each layer has named dimensions: `point`, `Ti` (time), and `depth` or `height`
where they apply. A monthly climatology gives one representative day per month, hour by hour.

```@example get_started
soil_surface = output.soil_temperature[depth = 1]
fig = Figure(size = (800, 300)) # hide
axes = [Axis(fig[1, i]; title = name, xlabel = "Hour of the representative day", ylabel = i == 1 ? "Soil surface (°C)" : "") # hide
        for (i, name) in enumerate(("Palm Springs", "San Jacinto Peak"))] # hide
linkyaxes!(axes...) # hide
for (i, ax) in enumerate(axes) # hide
    for (m, label) in ((1, "January"), (7, "July")) # hide
        lines!(ax, 0:23, celsius.(collect(soil_surface[point = i])[(m - 1) * 24 .+ (1:24)]); label) # hide
    end # hide
    i == 1 && axislegend(ax; position = :lt) # hide
end # hide
fig # hide
```

The weather is the same coarse climatology at both sites. The difference comes from the elevation in the DEM,
which corrects the air temperature, pressure and radiation for each point:

```@example get_started
terrain(cache).elevation
```

```@example get_started
markdown_table(["", "Palm Springs", "San Jacinto Peak"], [ # hide
    ("hottest soil surface (°C)", maximum(celsius.(soil_surface[point = 1])), maximum(celsius.(soil_surface[point = 2]))), # hide
    ("deepest snow (cm)", maximum(ustrip.(u"cm", output.snow_depth[point = 1])), maximum(ustrip.(u"cm", output.snow_depth[point = 2]))), # hide
]) # hide
```

## A map

The same model over a grid is a [`MicroRasterProblem`](@ref): an area and a grid to solve on. Here the grid is that
of CRU CL 2.0 itself, 10 arcminutes, over southern California:

```@example get_started
using Rasters.Extents: Extent
area = Extent(X = (-118.5, -115.5), Y = (32.8, 35.3))
grid_model = MicroMapModel(; micro_model, dem_source = CRUCL2, weather_source = CRUCL2,
    soil_moisture_source = CPCSoil, surface_albedo_source = 0.15, roughness_height_source = 0.004u"m",
    compute_terrain = false)
grid = MicroRasterProblem(; model = grid_model, area,
    template = load_template(CRUCL2, area),
    dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31),
    soil_profile = example_soil_profile(depths),
    init = (; soil_moisture = fill(0.2, length(depths)), snow_depth = 0.0u"cm"),
)
maps = solve(grid)
nothing # hide
```

```@example get_started
hottest = maximum(maps.soil_temperature[depth = 1]; dims = Ti)[Ti = 1]
snow = maximum(maps.snow_depth; dims = Ti)[Ti = 1]
map_figure(( # hide
    (celsius.(hottest), "Hottest soil surface (°C)", :thermal, "°C"), # hide
    (ustrip.(u"cm", snow), "Deepest snow (cm)", Reverse(:ice), "cm"), # hide
);  sites = ("Palm Springs" => sites[1],)) # hide
```

At 18 km a cell averages mountain and desert, so the cell that contains Palm Springs is far cooler, and snowier,
than the town. A finer DEM resolves them, at more cost, see [Maps](tutorials/maps.md).

## Another data set

Changing the weather is one keyword. WorldClim is another monthly climatology, for 1970–2000 where CRU CL 2.0 is
for 1961–1990, with different variables and grids:

```@example get_started
worldclim = MicroMapModel(; micro_model, dem_source = SRTM, weather_source = WorldClim{Climate},
    soil_moisture_source = CPCSoil, surface_albedo_source = 0.15, roughness_height_source = 0.004u"m")
other = solve(MicroVectorProblem(; model = worldclim, points = sites,
    dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31),
    soil_profile = example_soil_profile(depths),
    init = (; soil_moisture = fill(0.2, length(depths)), snow_depth = 0.0u"cm")))
markdown_table(["", "CRU CL 2.0", "WorldClim"], [ # hide
    ("hottest soil surface at Palm Springs (°C)", maximum(celsius.(soil_surface[point = 1])), maximum(celsius.(other.soil_temperature[point = 1, depth = 1]))), # hide
    ("deepest snow on San Jacinto Peak (cm)", maximum(ustrip.(u"cm", output.snow_depth[point = 2])), maximum(ustrip.(u"cm", other.snow_depth[point = 2]))), # hide
]) # hide
```

WorldClim publishes the elevation of its grid (`WorldClim{Elevation}` in RasterDataSources.jl), but the WorldClim
binding does not read it yet, so its air temperatures are not yet corrected to the elevation of each point, see
[Terrain and atmosphere](manual/terrain.md#Elevation). Daily and hourly data sets, such as
gridMET, NCEP and ERA5, are swapped in the same way and solve every day in turn, see
[Swapping data sets](tutorials/swap_data.md). See [One keyword per data set](manual/data_sources.md) for how this
works, and [Forcing data sets](manual/forcing_data.md) for the choices.

## Where next

- [Interactive app](interactive.md): pick a place on a map, swap data sets, plot the results and copy the script that
  reproduces them.
- [Points](tutorials/points.md) and [Maps](tutorials/maps.md) for more of each.
- [Built on Rasters and DimensionalData](manual/spatial_stack.md) for what can be done with the output.
- [Introduction](manual/introduction.md) for the design of the package.
