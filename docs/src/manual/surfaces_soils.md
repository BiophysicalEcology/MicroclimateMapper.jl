# Surfaces and soils

Besides weather and terrain, a microclimate needs the properties of the surface, its albedo and roughness, and
of the soil beneath it.

## Albedo and roughness

`surface_albedo_source` and `roughness_height_source` of the [`MicroMapModel`](@ref) each accept:

| Form | Example | Meaning |
|:--|:--|:--|
| a number | `0.15`, `0.004u"m"` | the same everywhere |
| a `Raster` | a map of albedo | resampled to the run's grid or points |
| a data-set type | a raster data set of the property | loaded and resampled |
| `nothing` | — | from `landcover_source`, by class, with the data set's default table |
| a NamedTuple | `(; Barren = 0.3, Snow_Ice = 0.8, …)` | from `landcover_source`, with your own values per class |

With land cover, each cell's value is the class values weighted by the fraction of each class in the cell
(`EarthEnv{LandCover}`, fractional) or taken from its single class (MODIS, categorical).

!!! warning "Placeholder values"
    The default albedo and roughness of each land-cover class are estimates, not yet taken from a published
    source. Supply your own NamedTuple where it matters.

## Vegetation

Under a canopy, the `shade` fraction can be replaced by Microclimate.jl's multilayer canopy model
(`MultilayerCanopy` in the `MicroModel`, see its [Canopy](https://biophysicalecology.github.io/Microclimate.jl/dev/manual/canopy)
page), which needs a description of the vegetation:

| Property | What it sets | Where it can come from |
|:--|:--|:--|
| canopy height (`canopy_height`) | the heights the canopy spans, and wind and air temperature within it | global canopy height maps, such as Lang et al. (2023) at 10 m; LiDAR, as the difference between surface and terrain models (see [Canopy on Saba](../tutorials/canopy.md)) |
| plant area index (PAI) and its profile with height (`plant_area_index`) | how much radiation each layer intercepts | leaf area index from satellite (MODIS, or Sentinel through the Copernicus WEkEO service), as an approximation: PAI also counts stems and branches |
| leaf angle distribution (χ, `canopy_projection_ratio` in `leaf_parameters`) | how steeply the leaves intercept the direct beam | values by vegetation type: near 0 for upright leaves (grasses, some conifers), 1 for leaves at all angles, larger for flat leaves |
| leaf reflectance and transmittance | how much shortwave is scattered within the canopy | per waveband; the model's defaults are 0.25 and 0.25 (`TwoStreamRadiation`) |
| woody area fraction | the part of PAI that is stems and branches, which intercept but do not transpire | by vegetation type |
| leaf size and stomatal conductance | leaf boundary layers and transpiration, hence leaf temperature (`leaf_parameters`, `stomatal_model`) | by vegetation type |

Ilya Maclean's R package [terravars](https://github.com/ilyamaclean/terravars) shows such layers in use: rasters
of PAI (from Sentinel leaf area index), canopy height (Lang et al. 2023) and leaf angle, from which it computes
two-stream canopy transmission (`twostream`), radiation through a three-dimensional canopy (`canopy3d`) and
the sky view through terrain and vegetation (`skyview`).

At present a run has one canopy, set on the point model, so a landscape of different vegetation is one run per
class of vegetation, see [Canopy on Saba](../tutorials/canopy.md). Reading these properties for each cell, from
rasters as for albedo and roughness, is planned, and not all of them are wired in yet.

## Soil profile

The `soil_profile` of a problem gives the bulk density, mineral density and hydraulic properties of each soil
layer, see Microclimate.jl's
[Soil thermal properties](https://biophysicalecology.github.io/Microclimate.jl/dev/manual/soil_thermal_properties)
and [Soil moisture](https://biophysicalecology.github.io/Microclimate.jl/dev/manual/soil_moisture).
`example_soil_profile(depths)` gives a typical loam. [`build_soil_profile`](@ref) builds one from soil texture
maps:

```julia
soil = build_soil_profile(SoilGrids, (-89.40, 43.07); depths)   # a point, anywhere
soil = build_soil_profile(SLGA, area; depths)                   # an area in Australia, averaged
problem = MicroVectorProblem(; model, points, dates, soil_profile = soil.soil_profile)
```

It reads clay, silt and sand fractions and bulk density at the depths of the data set,
[SoilGrids](https://soilgrids.org) (global, 250 m) or the
[Soil and Landscape Grid of Australia](https://www.clw.csiro.au/aclep/soilandlandscapegrid/) (90 m), interpolates
them onto the model's depths, and converts texture to the Campbell hydraulic parameters with a pedotransfer
function:

| `pedotransfer_model` | From |
|:--|:--|
| [`CosbyUnivariate`](@ref) | Cosby et al. (1984), table 5: from sand or clay fraction alone |
| [`CosbyMultivariate`](@ref) | Cosby et al. (1984), table 4: from sand, silt and clay |
| [`Campbell1985`](@ref) | Campbell (1985): from particle-size distribution and bulk density |

The field capacity and wilting point of each layer are returned for reference.

One soil profile is used for every point or cell of a run. An area is averaged into one profile. Soil that varies
from cell to cell is planned.

## Soil moisture

How soil moisture is treated is set by Microclimate.jl's `soil_moisture_strategy`:

- **prescribed** (the default): a time series of moisture is given, and only heat is solved;
- **dynamic**: infiltration, redistribution, evaporation and root uptake are solved from the rainfall, which
  costs more.

For prescribed moisture, MicroclimateMapper.jl takes the series from the weather data set where it has one
(`TerraClimate`), from `soil_moisture_source = CPCSoil` (a monthly climatology), or from `data = (; soil_moisture = …)`.

## Initial conditions

`init` on a problem gives the starting soil moisture, soil temperature and snow depth. Unset, soil temperature
starts at the mean air temperature of the first day. `init_source = ERA5` takes soil temperature and moisture from
ERA5 at the start date instead. Alternatively you can start the model earlier than your time of interest, i.e. inclue
a "spin-up" period, to allow time for the initial conditions' effects to fade (often six months to a year is needed).
See Microclimate.jl's [Initial conditions and burn-in](https://biophysicalecology.github.io/Microclimate.jl/dev/tutorials/initial_conditions_and_burnin).
