# Terrain and atmosphere

Weather data sets describe the air over cells kilometres wide. A point or a fine cell differs from that in two
ways: its elevation, which changes the air temperature and pressure, and its slope, aspect and horizon, which
change the radiation it receives. Both come from the `dem_source` of the [`MicroMapModel`](@ref).

```@setup terrain
using Main.FigureHelpers
```

## Slope, aspect and horizon

For every run the DEM is read over the area plus a margin (0.3° for SRTM, about 33 km, so that distant ridges
still count in the horizon), at its own resolution, and passed to
[Geomorphometry.jl](https://deltares.github.io/Geomorphometry.jl/stable/):

- **slope** and **aspect**, from the eight neighbours of each cell (Horn 1981), aspect clockwise from north;
- **horizon angles**, the elevation angle of the skyline in 32 (default) directions from north.

These are always computed on the DEM's own cells, the finest terrain available. Each cell of the run then gets a
`SolarTerrain` of [SolarRadiation.jl](https://biophysicalecology.github.io/SolarRadiation.jl/dev), from which
the direct and diffuse radiation on its surface are computed, with the sun hidden when it is below the horizon in
its direction.

### Grids

The run's grid is not necessarily the DEM's. A [`MicroRasterProblem`](@ref) solves on its `template`, which may be
any raster or a data-set type: the DEM itself, the DEM aggregated (`aggregate(mean, dem, 6)`), or the grid of the
weather (`load_template(CRUCL2, area)`). Everything is brought onto that template:

| Input | Onto the run's grid by |
|:--|:--|
| elevation, slope | the mean of the DEM cells within each run cell |
| aspect | a circular mean of the DEM cells (so that 350° and 10° average to north, not south) |
| horizon angles | the mean, direction by direction |
| weather | cubic-spline interpolation, so that it varies smoothly rather than in blocks of the weather grid |
| albedo, roughness, initial soil state | the mean |

With the DEM as template (or a finer DEM than the template) nothing is averaged: each cell has the DEM's own
elevation, slope, aspect and horizon. On a coarser template a cell's slope is the mean steepness of the ground
within it, not the slope of a smoothed surface, and its radiation is computed from that mean slope and aspect, an
approximation where the ground within the cell faces many ways.

A [`MicroVectorProblem`](@ref) has no template: each point takes the elevation, slope, aspect and horizon of the
DEM cell it falls in, and the weather of the weather cell it falls in.

SolarRadiation.jl's [Terrain](https://biophysicalecology.github.io/SolarRadiation.jl/dev/manual/terrain) page
explains the radiation, and its tutorial [Terrain: Saba](https://biophysicalecology.github.io/SolarRadiation.jl/dev/tutorials/saba)
does these steps by hand for an island in the Carribean: reversing the north-first rows, finding the cell size, rotating the
horizon directions to start at north. Here they are done for any DEM. [Terrain: Saba](../tutorials/saba.md) in
this documentation carries the Saba Island microclimate calculation through to soil temperatures.

The terrain of a run is kept on its cache, and can be passed to another run over the same grid to skip the DEM
and horizon calculations:

```julia
cache = init(problem)
again = MicroRasterProblem(; model = other_model, area, template, dates, soil_profile,
    data = (; terrain = terrain(cache)))
```

`compute_terrain = false` treats every cell as flat with an open horizon, which is quicker and suits analyses 
where flat ground calculations are sufficient.

## Elevation

The air temperature of the weather data is lapse-rate corrected from the elevation of its own grid to that of 
the point or cell of the DEM:

```@example terrain
using MicroclimateMapper, Unitful
markdown_table(["`lapse_rate_model`", "K per km", "At 1000 m above a 20 °C cell (°C)"],
    [("`$(nameof(typeof(lr)))`", ustrip(u"K/km", MicroclimateMapper.lapse_rate(lr)),
      ustrip(u"°C", MicroclimateMapper.lapse_adjust_temperature(lr, 20.0u"°C", 1000.0u"m")))
     for lr in (EnvironmentalLapseRate(), SaturatedAdiabaticLapseRate(), DryAdiabaticLapseRate(), CustomLapseRate(0.0077u"K/m"))])
```

The default is `EnvironmentalLapseRate()`. NicheMapR's `micro_global` used 7.7 K per km for maximum temperature
and 3.9 for minimum (`lapse_max`, `lapse_min`). Relative humidity is then recomputed at the corrected temperature,
the air pressure is set from the elevation, and wind speed is moved from the height of the data, usually 10 m, to
the model's reference height by a power law.

!!! warning "Only where the binding declares the grid's elevation"
    The correction needs the elevation of the weather grid: the surface the data set was made for. At present
    it is read only for `CRUCL2` and `BARRA`, whose bindings declare an elevation layer, so only their
    temperatures are corrected. For the others the temperature of the weather cell is used as it is, and only
    the pressure follows the DEM. Most of the others have such a surface: WorldClim publishes its elevation grid
    (`WorldClim{Elevation}`), as does gridMET (`GRIDMET{Elevation}`), the reanalyses have their model orography,
    and the gridded Australian and climatological data sets were interpolated over a DEM. Adding each to its
    binding is under way.

## Cloud and radiation

Data sets that give daily or monthly solar radiation but no cloud, such as `TerraClimate`, `WorldClim` and
`GRIDMET`, have their cloud cover estimated by comparing that radiation with the clear-sky radiation of
SolarRadiation.jl at the same place and day. The model then works from cloud cover, as Microclimate.jl expects.

With `solar_only = true` the microclimate is not solved at all, and the output is solar radiation alone, as
broadband or over wavebands (`SOLAR_PAR`, `SOLAR_UVB`, `SOLAR_NIR`), see [Outputs](outputs.md).
`cloud_correct_solar = true` reduces it by the cloud of the weather data; otherwise it is clear-sky.

## Aerosols

Clear-sky radiation depends on the aerosols in the air, see SolarRadiation.jl's
[Aerosols](https://biophysicalecology.github.io/SolarRadiation.jl/dev/manual/aerosols). The package includes a
reader for the Global Aerosol Data Set, `MicroclimateMapper.get_aerosol_optical_depth`, which gives the spectral
optical depth for a place, humidity and month. Using it for each cell of a run is planned; at present the
default atmosphere of SolarRadiation.jl is used everywhere.
