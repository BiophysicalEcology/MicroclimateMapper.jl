# Points

A [`MicroVectorProblem`](@ref) solves the microclimate at a list of points. For one point at a time, picked on a
map, the [interactive app](../interactive.md) does the same and writes the script. This tutorial runs three places that
appear elsewhere in these documentation sites: Madison, Wisconsin, the example site of Microclimate.jl; Palm
Springs, California, where the desert iguana of BiophysicalBehaviour.jl lives; and Mont Aigoual, the summit of the
Cévennes in southern France, which is also mapped in [Maps](maps.md).

```@setup points
using Main.FigureHelpers
using CairoMakie
```

```@example points
using MicroclimateMapper, Microclimate, RasterDataSources, Rasters, Dates, Unitful
using Microclimate: example_soil_profile, example_soil_properties_model, example_soil_hydraulic_model

depths = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 1.2]u"m"
micro_model = MicroModel(; depths, heights,
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
    snow_model = SnowModel())
model = MicroMapModel(; micro_model, dem_source = SRTM, weather_source = CRUCL2,
    soil_moisture_source = CPCSoil, surface_albedo_source = 0.15, roughness_height_source = 0.004u"m")
year = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)
start = (; soil_moisture = fill(0.2, length(depths)), snow_depth = 0.0u"cm")
nothing # hide
```

## Madison

Microclimate.jl's [Monthly workflow](https://biophysicalecology.github.io/Microclimate.jl/dev/tutorials/monthly_workflow)
drives Madison with `example_monthly_weather()`, the CRU CL 2.0 normals for the site, typed in. Here the same
normals are read from the data set:

```@example points
madison = solve(MicroVectorProblem(; model, points = [(-89.40123, 43.07305)], dates = year,
    soil_profile = example_soil_profile(depths), init = start))
nothing # hide
```

and the run of that tutorial, for comparison:

```@example points
reference = let model = MicroModel(; soil_properties_model = example_soil_properties_model(),
                                     soil_hydraulic_model = example_soil_hydraulic_model())
    solve(MicroProblem(model, MicroInputs(; site = example_site(), soil_profile = example_soil_profile(),
        environment_minmax = example_monthly_weather(), environment_daily = example_daily_environment(),
        environment_hourly = example_hourly_environment(), initial_soil_temperature = nothing,
        initial_soil_moisture = fill(0.42 * 0.25, length(model.depths)))))
end
fig = Figure(size = (800, 360)) # hide
for (k, (d_here, d_ref, label)) in enumerate(((1, 1, "surface"), (7, 13, "30 cm"))) # hide
    ax = Axis(fig[1, k]; title = label, xlabel = "Hour of the 12 representative days", ylabel = k == 1 ? "Soil temperature (°C)" : "") # hide
    lines!(ax, celsius.(reference.soil_temperature[:, d_ref]); label = "Microclimate.jl example", color = :grey60, linewidth = 3) # hide
    lines!(ax, celsius.(collect(madison.soil_temperature[point = 1, depth = d_here])); label = "MicroclimateMapper.jl, CRUCL2", color = :darkorange) # hide
end # hide
Legend(fig[2, 1:2], content(fig[1, 1]); orientation = :horizontal, framevisible = false) # hide
fig # hide
```

The two agree closely, as they should: the same climate, the same physics. They differ where the inputs do. Here
snow is modelled, and the flat winter traces at the surface are the snowpack, which holds the soil surface steady below freezing
while the bare soil of the example swings by 15 K a day; the elevation is that of SRTM at the point, and soil
moisture follows the CPC climatology, where the example holds it fixed.

## Selecting

The output is a `RasterStack`. Each layer is selected by named dimensions, not by column number:

```@example points
madison.soil_temperature[point = 1, depth = Near(0.1)]   # 10 cm (depths are in metres), every hour
nothing # hide
```

```@example points
summer = madison.air_temperature[point = 1, height = 2, Ti = Where(t -> month(t) == 7)]
extrema(summer)
```

## Palm Springs

The desert iguana of [A lizard's day](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/tutorials/lizard)
moves between the surface and a burrow. What it has to choose from, on the representative days of January and
July:

```@example points
palm_springs = solve(MicroVectorProblem(; model, points = [(-116.545, 33.830)], dates = year,
    soil_profile = example_soil_profile(depths), init = start))
fig = Figure(size = (800, 320)) # hide
for (k, (m, name)) in enumerate(((1, "January"), (7, "July"))) # hide
    ax = Axis(fig[1, k]; title = name, xlabel = "Hour", ylabel = k == 1 ? "Soil temperature (°C)" : "") # hide
    hours = (m - 1) * 24 .+ (1:24) # hide
    for (d, label) in ((1, "surface"), (4, "10 cm"), (7, "30 cm"), (9, "1 m")) # hide
        lines!(ax, 0:23, celsius.(collect(palm_springs.soil_temperature[point = 1, depth = d])[hours]); label) # hide
    end # hide
end # hide
Legend(fig[1, 3], content(fig[1, 1]); framevisible = false) # hide
fig # hide
```

BiophysicalBehaviour.jl takes microclimates like these through its `AvailableEnvironments`, see
[Activity and available environments](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/environments#From-gridded-climate-data).

## Mont Aigoual

Three sites on Mont Aigoual lie within 4 km of each other: a valley at Camprieu, a south-facing slope and the
summit. They share one cell of the weather data, so their differences come from elevation, slope and aspect:

```@example points
aigoual = solve(MicroVectorProblem(; model,
    points = [(3.520, 44.143), (3.582, 44.115), (3.581, 44.122)], dates = year,
    soil_profile = example_soil_profile(depths), init = start))
names = ["Valley", "South slope", "Summit"] # hide
fig = Figure(size = (800, 320)) # hide
ax = Axis(fig[1, 1]; xlabel = "Hour of the 12 representative days", ylabel = "Snow depth (cm)") # hide
for i in 1:3 # hide
    lines!(ax, ustrip.(u"cm", collect(aigoual.snow_depth[point = i])); label = names[i]) # hide
end # hide
Legend(fig[1, 3], ax; framevisible = false) # hide
ax2 = Axis(fig[1, 2]; xlabel = "Hour of the 12 representative days", ylabel = "Soil surface (°C)") # hide
for i in 1:3 # hide
    lines!(ax2, celsius.(collect(aigoual.soil_temperature[point = i, depth = 1])); label = names[i]) # hide
end # hide
fig # hide
```

A monthly climatology has no individual storms: each representative day carries the snow of the one before, so
the snowpack builds and melts smoothly with the months. Daily or hourly weather, such as the NCEP or ERA5 reanalyses,
gives real storms, see [Swapping data sets](swap_data.md).

## Points far apart

Points of one problem share one DEM, loaded over the box that contains them all. For points on different
continents, use a coarse DEM, such as `dem_source = CRUCL2` with `compute_terrain = false`, or one problem per
region, as here.
