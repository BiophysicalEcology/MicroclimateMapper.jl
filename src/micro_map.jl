# micro_map.jl
#
# Grid-mode microclimate driver. `MicroMapProblem` pairs a `MicroMapModel`
# with a rectangular spatial extent and a year range; `init` loads forcings
# at native resolution, resamples DEM-derived terrain (slope, aspect,
# horizon angles) onto the weather grid as the spatial template, and hands
# off to the shared scaffolding in `micro_common.jl`.
#
# `data::NamedTuple` on the problem provides two layers of override:
#   * collection-level — `weather::RasterStack`, `landcover::RasterStack|Raster`,
#     replacing the corresponding source loader output entirely
#   * single-layer — `dem`, `surface_albedo`, `roughness_height`, or any
#     canonical weather variable name (e.g. `vapour_pressure_deficit`,
#     `mean_temperature`), as a Raster in canonical units; resampled to the
#     run template at `init` time

"""
    MicroMapProblem(; model, area, years, init=nothing, data=(;))

Concrete grid-microclimate run: pairs a `MicroMapModel` with a spatial
extent, time range, initial conditions, and any per-run data overrides.

- `model::MicroMapModel`
- `area::Extents.Extent` — spatial extent (X = longitude, Y = latitude)
- `years::AbstractRange`
- `init` — initial conditions; default fills `soil_moisture` from
  `model.micro_model.depths` and lets `soil_temperature` fall back to the
  day-mean reference air temperature
- `data::NamedTuple` — overrides for any of:
    * collection sources: `weather` (RasterStack), `landcover` (Stack/Raster)
    * single layers: `dem`, `surface_albedo`, `roughness_height`, or any
      canonical weather-variable name (e.g. `vapour_pressure_deficit`,
      `mean_temperature`, `cloud_cover`). Each as a `Raster` in canonical
      units; resampled to the run template automatically.
"""
@kwdef struct MicroMapProblem{M<:MicroMapModel,A,Y,IT,D<:NamedTuple}
    model::M
    area::A
    years::Y
    init::IT = nothing
    data::D  = (;)
end

# ---------------------------------------------------------------------------
# Grid-mode spatial extractors
# ---------------------------------------------------------------------------

# Resample a native-resolution surface property onto the grid template.
# Scalar / Raster cases are handled here; the landcover-weighted path runs
# `_resolve_surface_native` first to compute a native Raster, then passes
# it through this same extractor.
_to_grid_template(value, template) =
    Raster(Fill(value, size(template, X), size(template, Y)),
           dims(template, (X, Y)))
_to_grid_template(r::Raster, template) = Rasters.resample(r; to = template)
@inline function _to_grid_template(r::AbstractArray{<:Quantity}, template)
    u = unit(eltype(r))
    stripped = ustrip.(u, r)
    resampled = Rasters.resample(rebuild(r, stripped); to = template)
    return rebuild(resampled, parent(resampled) .* u)
end

function _resolve_surface_grid(data_override, source, landcover_source, template, area, default_fn)
    native = data_override === nothing ?
        _resolve_surface_native(source, landcover_source, area, default_fn) :
        data_override
    return _to_grid_template(native, template)
end

# Pre-resample every canonical-override raster onto the run template once,
# so the per-pixel hot path is just an `Ti(k)` slice fetch.
function _resample_canonical_overrides(canonical_data::NamedTuple, template)
    return map(canonical_data) do raster
        _to_grid_template(raster, template)
    end
end

# ---------------------------------------------------------------------------
# CommonSolve.init
# ---------------------------------------------------------------------------

"""
    CommonSolve.init(problem::MicroMapProblem) -> MicroMapCache

Load all per-run data (weather, DEM, landcover, surface properties, canonical
overrides), pre-resolve every grid onto the run template, allocate the
worker-cache pool, and prepare for `solve!`.
"""
function CommonSolve.init(problem::MicroMapProblem)
    (; model, area, years, data) = problem
    (; dem_source, weather_source, landcover_source,
       surface_albedo_source, roughness_height_source) = model

    init_inputs = _resolve_init(problem.init, model.micro_model)
    soil_moisture_available = _has_canonical_input(:soil_moisture, weather_source, data)

    # Use the weather source as the spatial template; terrain stays at the
    # DEM's native resolution and is aggregated to the template inside
    # `compute_terrain_grids` so slope / aspect / horizon angles reflect
    # fine-scale relief rather than the smeared-out coarse-cell average.
    weather    = _resolve_weather(data, weather_source, area, years)
    dem_native = _resolve_dem(data, dem_source, area)
    template_2d = view(first(values(weather)), Ti(1))
    template    = first(values(weather))
    terrain     = compute_terrain_grids(dem_native;
        template = template_2d, n_horizon_angles = N_HORIZON_ANGLES)

    albedo_data    = get(data, :surface_albedo,   nothing)
    roughness_data = get(data, :roughness_height, nothing)
    albedo_grid    = _resolve_surface_grid(albedo_data, surface_albedo_source,
        landcover_source, template, area, default_landcover_albedo)
    roughness_grid = _resolve_surface_grid(roughness_data, roughness_height_source,
        landcover_source, template, area, default_landcover_roughness)

    canonical_overrides = _resample_canonical_overrides(_canonical_data(data), template)

    cloud_constants = _build_cloud_constants()
    build_inputs, cache_pool = _build_inputs_and_pool(;
        model, weather_source, weather, terrain,
        albedo_grid, roughness_grid, canonical_overrides,
        init_inputs, soil_moisture_available, years, cloud_constants,
    )

    return MicroMapCache(
        problem, weather, terrain, albedo_grid, roughness_grid,
        canonical_overrides, cache_pool,
        (; init_inputs, build_inputs),
        cloud_constants,
    )
end

"""
    CommonSolve.solve(problem::MicroMapProblem) -> RasterStack

Shortcut for `solve!(init(problem))`.
"""
CommonSolve.solve(problem::MicroMapProblem) = solve!(init(problem))

# ---------------------------------------------------------------------------
# Show methods
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", problem::MicroMapProblem)
    println(io, "MicroMapProblem")
    println(io, "  area:   ", problem.area)
    println(io, "  years:  ", problem.years)
    println(io)
    println(io, "forcings:")
    _print_role_table(io, _role_statuses(problem))
    _print_overrides(io, problem.data, "")
    _print_problem_init(io, problem)
end

function Base.show(io::IO, ::MIME"text/plain", cache::MicroMapCache{<:MicroMapProblem})
    problem = cache.problem
    nx, ny = size(cache.terrain.elevation)
    println(io, "MicroMapCache  (", nx, "×", ny, " grid, ", nx * ny, " pixels)")
    println(io, "  cache_pool:  ", cache.cache_pool.sz_max, " workers")
    println(io)
    println(io, "forcings:")
    _print_role_table(io, _role_statuses(problem))
    _print_overrides(io, problem.data, "")
end
