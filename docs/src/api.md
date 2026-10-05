# API

```@meta
CurrentModule = MicroclimateMapper
```

## Model and problems

```@docs
MicroMapModel
MicroVectorProblem
MicroRasterProblem
MicroMapCache
```

## Solving

```@docs
CommonSolve.init(::MicroRasterProblem)
CommonSolve.init(::MicroVectorProblem)
CommonSolve.solve(::MicroRasterProblem)
CommonSolve.solve(::MicroVectorProblem)
CommonSolve.solve!(::MicroMapCache)
terrain
prefetch_weather!
```

## Terrain and lapse rates

```@docs
load_template
LapseRate
EnvironmentalLapseRate
DryAdiabaticLapseRate
SaturatedAdiabaticLapseRate
CustomLapseRate
```

## Soils

```@docs
build_soil_profile
PedotransferModel
CosbyUnivariate
CosbyMultivariate
Campbell1985
```

## Outputs

```@docs
LayerSpec
canonical_unit
strip_to_canonical
SolarOutputLayer
SOLAR_BROADBAND
SOLAR_PAR
SOLAR_UVB
SOLAR_NIR
```

## Places

```@docs
geocode
GeocodeResult
```

## Data-set loaders

```@docs
MicroclimateMapper.PointQuery
```

## Interactive app

See [Interactive](interactive.md).

```@docs
app
```
