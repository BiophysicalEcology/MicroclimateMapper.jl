```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "MicroclimateMapper.jl"
  text: "Microclimates anywhere, from any data"
  tagline: "Microclimate.jl over points and grids, driven by global and regional weather, terrain, soil and land-cover data sets, chosen with one keyword."
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/BiophysicalEcology/MicroclimateMapper.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - title: 🔑 One keyword per data set
    details: CRU CL 2.0, WorldClim, CHELSA, TerraClimate, gridMET, NCEP, ERA5, AWAP, SILO and BARRA, swapped with <a class="highlight-link">weather_source</a> and nothing else.
    link: /manual/data_sources
  - title: 🗺️ Points or maps
    details: The same model for one point, many points or every cell of a grid, returned as a RasterStack with named dimensions.
    link: /tutorials/maps
  - title: ⛰️ Terrain from any DEM
    details: Slope, aspect and horizon angles from <a class="highlight-link">Geomorphometry.jl</a>, and weather corrected for elevation.
    link: /manual/terrain
  - title: 🧩 Built on Rasters.jl
    details: Data from <a class="highlight-link">RasterDataSources.jl</a>, held and analysed with Rasters.jl and DimensionalData.jl.
    link: /manual/spatial_stack
  - title: 🌱 Soils and land cover
    details: Soil profiles from SoilGrids and SLGA texture, albedo and roughness from land-cover classes.
    link: /manual/surfaces_soils
  - title: 🌳 Canopies
    details: Microclimate.jl's multilayer canopy driven across a landscape.
    link: /tutorials/canopy
  - title: 🦎 From NicheMapR
    details: The Julia counterpart of micro_global, micro_usa, micro_ncep and their relatives.
    link: /manual/nichemapr
---
```

## How to install MicroclimateMapper.jl?

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/BiophysicalEcology/MicroclimateMapper.jl")
```

## Manual

MicroclimateMapper.jl runs [Microclimate.jl](https://biophysicalecology.github.io/Microclimate.jl/dev) over
points and grids, taking its weather, terrain, soil and surface inputs from spatial data sets. See the
[Introduction](manual/introduction.md) for the design, and [For NicheMapR users](manual/nichemapr.md) for its
origins in [NicheMapR](https://github.com/mrke/NicheMapR).

It is part of the [BiophysicalEcology](https://github.com/BiophysicalEcology) ecosystem for mechanistic niche
modelling, and builds on [RasterDataSources.jl](https://ecojulia.github.io/RasterDataSources.jl/stable),
[Rasters.jl](https://rafaqz.github.io/Rasters.jl/stable),
[DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl/stable/) and
[Geomorphometry.jl](https://deltares.github.io/Geomorphometry.jl/stable/).
