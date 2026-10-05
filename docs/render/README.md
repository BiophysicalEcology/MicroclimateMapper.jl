# Render scripts

Some figures in the documentation come from runs too slow, or data too large, for the documentation build: grids
of thousands of cells, and data sets that are large downloads or need credentials. These scripts make them. Each
writes its figures and tables into `docs/src/assets/generated/`, which are committed, and the pages show the code
that produced them.

| Script | Pages | Data | Time here |
|:--|:--|:--|:--|
| `swap_data.jl [Madison] [Melbourne]` | Swapping data sets | CRUCL2, WorldClim, TerraClimate 2000, SILO 2025; one untimed warm-up run first | ~5 min |
| `maps.jl [pssj] [aigoual]` | Maps | SRTM, CRUCL2 | ~30 min |
| `saba.jl` | Terrain: Saba | the Saba LiDAR DEM, TerraClimate 2000 | ~35 min |
| `canopy.jl` | Canopy on Saba | the Saba DEM and DSM, TerraClimate 2000; soil, leaf and canopy air temperatures | ~3 min |
| `soils.jl` | Soils and land cover | SLGA; SILO 2025 with dynamic soil moisture (`soils_moisture.png`, not yet used: see the script) | ~10 min |
| `nichemapr.jl` | For NicheMapR users | TerraClimate 2000, `test/data/micro_terra` | ~2 min |

Run them with the documentation's environment and all threads, from the repository root:

```
julia --project=docs --threads=auto docs/render/maps.jl
```

Data are read from `ENV["RASTERDATASOURCES_PATH"]`, `c:/Spatial_Data/` if unset. `maps.jl` keeps the Mont Aigoual
output in `docs/render/cache/` (not committed) so its figures can be redrawn without solving again. The times are
from a 16-thread laptop in October 2026.
