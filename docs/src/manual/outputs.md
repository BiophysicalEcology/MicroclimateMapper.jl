# Outputs

`solve` returns a `RasterStack`: one layer per output variable, all on the same spatial and time dimensions.

| Dimension | In | Values |
|:--|:--|:--|
| `point` | a [`MicroVectorProblem`](@ref) | the `(longitude, latitude)` of each point |
| `X`, `Y` | a [`MicroRasterProblem`](@ref) | the run's grid |
| `Ti` | all | the hours solved: every hour of every day for daily data, every hour of 12 representative days for monthly |
| `depth` | soil layers | the `depths` of the `MicroModel` |
| `height` | air layers | the `heights` of the `MicroModel` |

Cells outside an `area` polygon, or with no data, are `missing`.

## Layers

Which layers are written is set by `output_layers` of the [`MicroMapModel`](@ref), a tuple of
[`LayerSpec`](@ref)s. Each names a field of Microclimate.jl's result and says what kind it is:

| Kind | Shape per point or cell | Default layers |
|:--|:--|:--|
| `:soil` | time × depth | `soil_temperature`, `soil_moisture` |
| `:profile` | time × height | `air_temperature`, `relative_humidity`, `wind_speed` |
| `:scalar` | time | `global_radiation`, `sky_temperature`, `snow_depth`, `ground_surface_water`, `ground_dew`, `ground_frost`, `ground_standing_dew`, `ground_standing_frost` |
| `:canopy` | time × canopy layer | none by default; with a multilayer canopy, e.g. `LayerSpec(:leaf_temperature, :canopy)` |

Writing fewer layers saves memory, which limits grid size before time does:

```julia
model = MicroMapModel(; micro_model, dem_source = SRTM, weather_source = CRUCL2,
    output_layers = (LayerSpec(:soil_temperature, :soil), LayerSpec(:snow_depth, :scalar)))
```

Where two kinds share a field name, a third argument names the field and the first names the layer:
`LayerSpec(:canopy_air_temperature, :canopy, :air_temperature)`.

## Units

Every layer carries units. For writing to NetCDF or another format without units, [`strip_to_canonical`](@ref)
converts each layer to the unit given by [`canonical_unit`](@ref) and strips it:

```julia
write("run.nc", strip_to_canonical(output))
```

| Layers | Canonical unit |
|:--|:--|
| temperatures | °C |
| `soil_moisture` | m³ m⁻³ |
| `relative_humidity` | % |
| `wind_speed` | m s⁻¹ |
| `global_radiation` | W m⁻² |
| `snow_depth` | cm |
| water at the surface, dew and frost | kg m⁻² |

## Solar radiation

`solar_output_layers` adds layers of solar radiation, alongside the microclimate or, with `solar_only = true`,
instead of it. Each is a [`SolarOutputLayer`](@ref): a component (global, direct or diffuse on a horizontal
surface, or global on the sloping surface of the cell) and optionally a waveband. Four are predefined:

| Constant | Layer | Waveband |
|:--|:--|:--|
| [`SOLAR_BROADBAND`](@ref) | `global_radiation` | all |
| [`SOLAR_PAR`](@ref) | `par` | 400–700 nm |
| [`SOLAR_UVB`](@ref) | `uv_b` | 290–315 nm |
| [`SOLAR_NIR`](@ref) | `nir` | 700–4000 nm |

Wavebands come from the spectral model of SolarRadiation.jl, see its
[Atmosphere](https://biophysicalecology.github.io/SolarRadiation.jl/dev/manual/atmosphere) page.
