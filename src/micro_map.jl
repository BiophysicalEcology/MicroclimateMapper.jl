# micro_map.jl
#
# Grid-mode microclimate driver. `MicroMapProblem` pairs a `MicroMapModel`
# with a rectangular spatial extent, a year range, and an explicit
# `template` grid (a Raster or an RDS source loaded over `area`); `init`
# loads forcings at native resolution, resamples weather and DEM-derived
# terrain (slope, aspect, horizon angles) onto the template, and hands
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
    MicroMapProblem(; model, area, years, template, init=nothing, data=(;))

Concrete grid-microclimate run: pairs a `MicroMapModel` with a spatial
extent, time range, initial conditions, and any per-run data overrides.

- `model::MicroMapModel`
- `area::Extents.Extent` — spatial extent (X = longitude, Y = latitude)
- `years::AbstractRange`
- `template` — spatial grid the run executes on. Required; no fallback.
    * `Type{<:RasterDataSource}` (e.g. `SRTM`) — load that dataset over
      `area` and use it as the grid; weather is resampled onto it.
    * `Raster` — use the supplied raster's `X`/`Y` lookup as the grid;
      weather is resampled onto it. Only `X`/`Y` dims are read.
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
@kwdef struct MicroMapProblem{M<:MicroMapModel,A,Y,IT,D<:NamedTuple,T}
    model::M
    area::A
    years::Y
    template::T
    init::IT    = nothing
    data::D     = (;)
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

# Resolve the user-supplied `MicroMapProblem.template` into a 2-D Raster.
#   Raster              → use directly (slice off `Ti` if it carries one).
#   Type{<:RDS source}  → load over `area` and use that as the template.
_resolve_template(r::Raster, _area) = hasdim(r, Ti) ? view(r, Ti(1)) : r
# Only the X/Y dims are read downstream (via `Rasters.resample(; to=template)`
# and `dims(template, (X, Y))`), so keep the raster lazy — never materialise.
# Pass `extent` only if the source advertises it in `getraster_keywords`
# (e.g. SRTM mosaics tile selection on it); otherwise load the whole
# dataset lazily and let `crop` trim down to `area`.
function _resolve_template(source::Type{<:RasterDataSources.RasterDataSource}, area::Extent)
    kwargs = :extent in RasterDataSources.getraster_keywords(source) ?
        (; extent = area, lazy = true) : (; lazy = true)
    return crop(Raster(source; kwargs...); to = area, touches = true)
end

# One-line description of the template for `show`.
_template_label(r::Raster) = "user Raster ($(size(r, X))×$(size(r, Y)))"
_template_label(source::Type{<:RasterDataSources.RasterDataSource}) = string(source)

# Degree buffer applied to the weather load area. Needs to be wide enough
# that the cubic-spline resample kernel (4-cell footprint) has source
# cells beyond the requested template extent.
#
# The default (0.5°) covers 1/24° (TerraClimate), 1/120° (CHELSA,
# WorldClim), and 0.25° (ERA5) native grids. Coarser sources override
# this — NCEP's T62 Gaussian grid is ~1.9°, so 4° gives the spline two
# source cells beyond the template extent in every direction.
weather_area_buffer(::Type{<:RasterDataSources.RasterDataSource}) = 0.5
weather_area_buffer(::Type{<:NCEP}) = 4.0

# Resample every layer of a 3-D `(X, Y, Ti)` weather stack onto a 2-D
# `(X, Y)` template — used when the run grid is finer than the native
# weather grid (e.g. SRTM-resolution template over TerraClimate weather).
# `:cubicspline` interpolates continuously across native cell boundaries
# so the upsampled weather varies smoothly with terrain rather than
# stepping in blocks aligned to the coarse weather grid (GDAL's default
# is nearest, which would produce visible block artifacts).
function _resample_weather_to_template(weather::RasterStack, template_2d)
    target_crs = crs(template_2d)
    names = keys(weather)
    resampled = map(names) do n
        layer = _regularise_for_resample(getproperty(weather, n), target_crs)
        Rasters.resample(layer; to = template_2d, method = :cubicspline)
    end
    return RasterStack(NamedTuple{names}(resampled))
end

# GDAL's resample backend needs a `Projected` lookup it can warp from, but
# Gaussian-grid sources (NCEP) come in as `Mapped` with `Irregular` Y
# bounds, which `convertlookup(::Projected, ::Mapped{Irregular})` can't
# materialise. Locally — over the buffered template extent — Gaussian
# spacing is uniform enough that rebuilding Y as a `Regular`-span
# `Sampled` lookup at the mean step gives identical interpolation
# results. We also rebuild X as plain `Sampled`+target CRS in the same
# pass so neither axis goes through the `Mapped`→`Projected` reproject
# path (which would otherwise require loading Proj).
function _regularise_for_resample(layer::AbstractRaster, target_crs)
    _is_irregular(lookup(layer, Y)) || return layer
    new_x = _regular_sampled(lookup(layer, X))
    new_y = _regular_sampled(lookup(layer, Y))
    return setcrs(Rasters.set(layer, X => new_x, Y => new_y), target_crs)
end

@inline _is_irregular(lk) = span(lk) isa Irregular

# Build a `Sampled{Regular}` lookup whose values land at the mean spacing
# between the existing coordinates. Reversed inputs (descending Y) keep
# their order. `Intervals(Center)` matches GDAL's cell-centre convention.
function _regular_sampled(lk)
    vals = collect(lk)
    n = length(vals)
    if n < 2
        return Sampled(Float64.(vals); sampling = Intervals(Center()),
            order = ForwardOrdered(), span = Regular(1.0))
    end
    step = (Float64(vals[end]) - Float64(vals[1])) / (n - 1)
    rng  = range(Float64(vals[1]); step = step, length = n)
    ord  = step < 0 ? ReverseOrdered() : ForwardOrdered()
    return Sampled(rng; sampling = Intervals(Center()),
        order = ord, span = Regular(step))
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

    # Resolve the spatial template (Raster or RDS source loaded over
    # `area`), then resample the native-resolution weather onto it so
    # per-pixel `weather[I..., Ti(k)]` indexing matches the run grid.
    # Terrain stays at the DEM's native resolution and is aggregated to
    # the template inside `compute_terrain_grids`, so slope / aspect /
    # horizon angles reflect fine-scale relief regardless of template.
    # Weather is loaded over a buffered area so the cubic-spline kernel
    # used by `_resample_weather_to_template` has source cells beyond
    # the template edges; without the buffer the spline degenerates to
    # nearest at the boundary.
    buffer_deg = weather_area_buffer(weather_source)
    weather_area = Extents.buffer(area,
        (X = buffer_deg, Y = buffer_deg))
    weather     = _resolve_weather(data, weather_source, weather_area, years)
    dem_native  = _resolve_dem(data, dem_source, area)
    template_2d = _resolve_template(problem.template, area)
    weather     = _resample_weather_to_template(weather, template_2d)
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
    println(io, "  area:     ", problem.area)
    println(io, "  years:    ", problem.years)
    println(io, "  template: ", _template_label(problem.template))
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
