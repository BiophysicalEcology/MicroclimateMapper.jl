# micro_map.jl
#
# Declarative grid-microclimate API:
#
#   MicroMapModel   — declarative model description (sources, output layout,
#                     lapse-rate model, inner per-pixel `MicroModel`)
#   MicroMapProblem — concrete run (model + spatial extent + years + initial
#                     conditions + per-run data overrides)
#   MicroMapCache   — workspace, built by `init(problem)` and consumed by
#                     `solve!(cache)`
#
# The interface mirrors Microclimate.jl: `solve(problem)` does everything;
# `init(problem)` then `solve!(cache)` gives access to the workspace.
#
# `data::NamedTuple` on the problem provides two layers of override:
#   * collection-level — `weather::RasterStack`, `landcover::RasterStack|Raster`,
#     replacing the corresponding source loader output entirely
#   * single-layer — `dem`, `surface_albedo`, `roughness_height`, or any
#     canonical weather variable name (e.g. `vapour_pressure_deficit`,
#     `mean_temperature`), as a Raster in canonical units; resampled to the
#     run template at `init` time

# ---------------------------------------------------------------------------
# Output layer specs and per-pixel write
# ---------------------------------------------------------------------------

# A `LayerSpec{Name, Kind}` encodes both the result-field name and its storage
# kind in its type, so each unrolled iteration specialises on a concrete type.
# `Kind` is one of:
#   :soil    — `(X, Y, Ti, Depth)`,  sourced from `result.<Name>`
#   :profile — `(X, Y, Ti, Height)`, sourced from `result.profile.<Name>`
#   :scalar  — `(X, Y, Ti)`,         sourced from `result.<Name>`
# Number of azimuth directions used for horizon-angle computation. A
# power of two makes the cardinal directions land on sample points.
const N_HORIZON_ANGLES = 32

struct LayerSpec{Name, Kind} end
LayerSpec(name::Symbol, kind::Symbol) = LayerSpec{name, kind}()

@inline _layer_name(::LayerSpec{N}) where N = N
@inline _layer_source(result, ::LayerSpec{N, :profile}) where N = getproperty(result.profile, N)
@inline _layer_source(result, ::LayerSpec{N}) where N = getproperty(result, N)

const _DEFAULT_OUTPUT_LAYERS = (
    LayerSpec(:soil_temperature,  :soil),
    LayerSpec(:soil_moisture,     :soil),
    LayerSpec(:air_temperature,   :profile),
    LayerSpec(:relative_humidity, :profile),
    LayerSpec(:wind_speed,        :profile),
    LayerSpec(:surface_water,     :scalar),
    LayerSpec(:global_radiation,  :scalar),
    LayerSpec(:sky_temperature,   :scalar),
)

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
    MicroMapModel(; micro_model, dem_source, weather_source,
                    landcover_source=nothing,
                    surface_albedo_source=nothing,
                    roughness_height_source=nothing,
                    output_layers=_DEFAULT_OUTPUT_LAYERS,
                    lapse_rate_model=EnvironmentalLapseRate())

Declarative description of a grid-microclimate run. Constant across spatial
extents and time ranges — pair with a `MicroMapProblem` to actually run.

- `micro_model::MicroModel` — inner per-pixel physics model
- `dem_source` — DEM data-source type (e.g. `SRTM`)
- `weather_source` — weather data-source type (e.g. `TerraClimate{Historical}`)
- `landcover_source` — optional land-cover data-source type
  (e.g. `EarthEnv{LandCover}`)
- `surface_albedo_source`, `roughness_height_source` — accepted forms:
    * `nothing` (use `landcover_source`'s default class→value table)
    * a class→value `NamedTuple` (overrides the dataset default)
    * a scalar (broadcast)
    * a `Raster` (resampled to the weather template)
    * a `Type{<:RasterDataSource}` (load a property raster directly)
- `output_layers` — tuple of `LayerSpec`s controlling which output rasters
  are materialised
- `lapse_rate_model::LapseRate` — atmospheric lapse-rate model for
  elevation-correcting weather data
"""
@kwdef struct MicroMapModel{MM,DS,WS,LCS,SAS,RHS,OL,LRT}
    micro_model::MM
    dem_source::DS
    weather_source::WS
    landcover_source::LCS        = nothing
    surface_albedo_source::SAS   = nothing
    roughness_height_source::RHS = nothing
    output_layers::OL            = _DEFAULT_OUTPUT_LAYERS
    lapse_rate_model::LRT         = EnvironmentalLapseRate()
end

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

# Per-run workspace built by `init(problem)`. Holds the loaded forcings,
# pre-resolved surface grids, and a Channel pool of per-worker inner
# `MicroCache`s + scratch buffers. The output `RasterStack` is allocated
# fresh on every `solve!` call — it does not live on the cache.
mutable struct MicroMapCache{P<:MicroMapProblem,W,T,A,R,CO,POOL,SC,CC}
    problem::P
    weather::W                       # RasterStack
    terrain::T                       # RasterStack
    albedo_grid::A                   # Raster (X, Y)
    roughness_grid::R                # Raster (X, Y)
    canonical_overrides::CO          # NamedTuple of pre-resampled Rasters
    cache_pool::POOL                 # Channel{(; micro::MicroCache, scratch)}
    init_inputs::SC                  # Resolved initial conditions
    cloud_constants::CC              # Shared SolarRadiation constants
end

# ---------------------------------------------------------------------------
# Data override resolution
# ---------------------------------------------------------------------------

# Keys in `problem.data` that target collection sources or named single
# layers; anything else is interpreted as a canonical weather variable.
const _SPECIAL_DATA_KEYS = (
    :weather, :landcover, :dem, :surface_albedo, :roughness_height,
)

@inline _is_special_key(k::Symbol) = k in _SPECIAL_DATA_KEYS

# Filter `data` down to the canonical-weather-variable overrides
# (anything not in `_SPECIAL_DATA_KEYS`).
function _canonical_data(data::NamedTuple)
    ks = filter(k -> !_is_special_key(k), keys(data))
    return NamedTuple{ks}(map(k -> data[k], ks))
end

# Pre-resample every canonical-override raster onto the run template once,
# so the per-pixel hot path is just an `Ti(k)` slice fetch.
function _resample_canonical_overrides(canonical_data::NamedTuple, template)
    return map(canonical_data) do raster
        _resample_canonical_override(raster, template)
    end
end

@inline function _resample_canonical_override(raster::AbstractArray{<:Quantity},
                                              template)
    # GDAL can't handle unitful arrays — strip, resample, re-attach.
    u = unit(eltype(raster))
    stripped = ustrip.(u, raster)
    resampled = Rasters.resample(rebuild(raster, stripped); to = template)
    return rebuild(resampled, parent(resampled) .* u)
end
@inline _resample_canonical_override(raster, template) =
    Rasters.resample(raster; to = template)

# Resolve initial conditions. Both keys default to `nothing`, meaning
# "fall back to the weather source / data override" — never fabricated.
_resolve_init(it::NamedTuple, _) = it
_resolve_init(::Nothing, _) =
    (soil_temperature = nothing, soil_moisture = nothing)

# Does the run actually have a `:soil_moisture` value to read? True iff
# the weather source declares it as a canonical variable or the user
# supplied `data.soil_moisture`.
function _has_canonical_input(name::Symbol, weather_source, data::NamedTuple)
    haskey(data, name) && return true
    weather_source === nothing && return false
    for var in weather_variables(weather_source)
        canonical_name(var) === name && return true
    end
    return false
end

# Per-pixel `initial_soil_moisture` resolution.
#  - User-supplied vector wins.
#  - Otherwise read the first timestep of `buffers.soil_moisture` and
#    broadcast across depths (TerraClimate-style rooting-zone average).
#  - If neither is available, refuse to fabricate.
_initial_soil_moisture(user::AbstractVector, _, _, _) = user
function _initial_soil_moisture(::Nothing, buffer, available::Bool, depths)
    available || error(
        "No soil_moisture available. Either: pass `init = (; soil_moisture = vec)`, " *
        "use a weather source that provides it as a canonical variable (e.g. TerraClimate), " *
        "or supply `data.soil_moisture` on the MicroMapProblem."
    )
    return fill(buffer[1], length(depths))
end

# ---------------------------------------------------------------------------
# Source resolution
# ---------------------------------------------------------------------------

# DEM: override beats source.
_resolve_dem(data::NamedTuple, source, area) =
    haskey(data, :dem) ? data.dem : _load_dem(source, area)

# Weather stack: override beats source.
_resolve_weather(data::NamedTuple, source, area, years) =
    haskey(data, :weather) ? data.weather : _load_weather(source, area, years)

# Surface-property resolution extended with two extra inputs:
#   * a `data` Raster override for that exact property (highest priority)
#   * a `Type{<:RasterDataSource}` as the model's `*_source` (load via RDS,
#     resample to template — no landcover weighting)
function _resolve_surface_grid(
    data_override, source, landcover_source, template, area, default_fn,
)
    data_override === nothing || return Rasters.resample(data_override; to = template)
    return _resolve_surface_property(source, landcover_source, template, area, default_fn)
end

# ---------------------------------------------------------------------------
# CommonSolve.init
# ---------------------------------------------------------------------------

"""
    CommonSolve.init(problem::MicroMapProblem) -> MicroMapCache

Load all per-run data (weather, DEM, landcover, surface properties, canonical
overrides), pre-resolve every grid onto the run template, allocate the
output stack, build a per-worker `MicroCache` pool, and solve the first
pixel (needed to size the output stack from its result eltypes).
"""
function CommonSolve.init(problem::MicroMapProblem)
    (; model, area, years, data) = problem
    (; micro_model, dem_source, weather_source, landcover_source,
       surface_albedo_source, roughness_height_source,
       output_layers, lapse_rate_model) = model

    init_inputs = _resolve_init(problem.init, micro_model)
    soil_moisture_available = _has_canonical_input(:soil_moisture, weather_source, data)

    # Load forcings (or take override). Use the weather source as the
    # spatial template; terrain stays at the DEM's native resolution and is
    # aggregated to the template inside `compute_terrain_grids` so slope /
    # aspect / horizon angles reflect fine-scale relief rather than the
    # smeared-out coarse-cell average.
    weather    = _resolve_weather(data, weather_source, area, years)
    dem_native = _resolve_dem(data, dem_source, area)
    template_2d = view(first(values(weather)), Ti(1))
    template   = first(values(weather))
    terrain    = compute_terrain_grids(dem_native;
        template = template_2d, n_horizon_angles = N_HORIZON_ANGLES)

    albedo_data    = get(data, :surface_albedo,   nothing)
    roughness_data = get(data, :roughness_height, nothing)
    albedo_grid    = _resolve_surface_grid(albedo_data, surface_albedo_source,
        landcover_source, template, area, default_landcover_albedo)
    roughness_grid = _resolve_surface_grid(roughness_data, roughness_height_source,
        landcover_source, template, area, default_landcover_roughness)

    canonical_overrides = _resample_canonical_overrides(_canonical_data(data), template)

    # Shared, immutable, build-once-for-the-run constants used by
    # `assemble_weather!`. Per-pixel elevation/pressure/lat/lon are swapped
    # into `flat_terrain_template` via `setproperties` inside the call.
    cloud_constants = (;
        solar_model = SolarProblem(),
        hours = collect(0.0:1.0:23.0),
        flat_terrain_template = SolarTerrain(;
            elevation = 0.0u"m", slope = 0.0u"°", aspect = 0.0u"°",
            horizon_angles = Fill(0.0u"°", N_HORIZON_ANGLES),
            albedo = 0.15,
            atmospheric_pressure = atmospheric_pressure(0.0u"m"),
            latitude = 0.0u"°", longitude = 0.0u"°",
        ),
    )

    # Per-pixel solar scratch sizing — `ndays` counts the main timesteps
    # the source provides (12 monthly representative days for monthly
    # sources, 365 calendar days for daily sources, etc.).
    ndays  = steps_per_year(temporal_resolution(weather_source)) * length(years)
    nsteps = length(cloud_constants.hours) * ndays
    nmax   = cloud_constants.solar_model.wavelength_count
    allocate_scratch() = (;
        weather = allocate_weather_buffers(weather_source, length(years)),
        solar = (;
            out = allocate_output_arrays(nsteps, ndays, nmax),
            buffers = allocate_buffers(nmax, cloud_constants.solar_model.diffuse_model),
        ),
        cloud_constants,
    )

    vapour_pressure_method = micro_model.vapour_pressure_equation

    function build_inputs(scratch, I::Tuple)
        horizon_angles = terrain.horizon_angles[I...]
        site = Site(;
            elevation = terrain.elevation[I...],
            slope = terrain.slope[I...],
            aspect = terrain.aspect[I...],
            latitude = terrain.latitude[I...],
            longitude = terrain.longitude[I...],
            atmospheric_pressure = terrain.atmospheric_pressure[I...],
            horizon_angles,
            sky_view_fraction = _sky_view_from_horizon(horizon_angles),
            albedo = albedo_grid[I...],
            roughness_height = roughness_grid[I...],
        )
        env = assemble_weather!(scratch, weather, weather_source, site, I;
            vapour_pressure_method, lapse_rate_model, canonical_overrides)
        initial_soil_moisture = _initial_soil_moisture(
            init_inputs.soil_moisture, scratch.weather.soil_moisture,
            soil_moisture_available, micro_model.depths,
        )
        MicroInputs(; site,
            env.environment_minmax, env.environment_daily, env.environment_hourly,
            initial_soil_temperature = init_inputs.soil_temperature,
            initial_soil_moisture,
        )
    end

    npixels = length(terrain.elevation)
    first_I = first(DimIndices(terrain.elevation))
    build_cache() = let scratch = allocate_scratch()
        (micro = CommonSolve.init(MicroProblem(micro_model, build_inputs(scratch, first_I))),
         scratch)
    end

    # Per-worker cache pool. The model is shared by reference (immutable);
    # each worker gets its own scratch and inner MicroCache. `init` only
    # allocates — no pixel is solved here.
    nworkers = min(Threads.nthreads(), npixels)
    proto = build_cache()
    cache_pool = Channel{typeof(proto)}(nworkers)
    put!(cache_pool, proto)
    for _ in 2:nworkers
        put!(cache_pool, build_cache())
    end

    return MicroMapCache(
        problem, weather, terrain, albedo_grid, roughness_grid,
        canonical_overrides, cache_pool,
        (; init_inputs, build_inputs),
        cloud_constants,
    )
end

# ---------------------------------------------------------------------------
# CommonSolve.solve! / solve
# ---------------------------------------------------------------------------

"""
    CommonSolve.solve!(cache::MicroMapCache) -> RasterStack
    CommonSolve.solve!(output::RasterStack, cache::MicroMapCache) -> RasterStack

Run the per-pixel microclimate solve across every grid pixel. Pixels are
dispatched over the cache pool with one `Threads.@spawn` worker per pool
slot.

The one-arg form allocates a fresh `RasterStack` to write into — sized
from the first pixel's result so the output's element types are concrete.
The two-arg form writes into the supplied `output` and returns it, for
callers that want to reuse storage across runs.
"""
function CommonSolve.solve!(cache::MicroMapCache)
    proto, first_I, first_result = _solve_proto_pixel!(cache)
    try
        layers = cache.problem.model.output_layers
        output = _allocate_output(cache.problem.model.micro_model,
            cache.terrain, first_result, layers, first(cache.problem.years))
        _write_output!(output, first_result, layers, first_I)
        return _solve_remaining!(output, cache, proto)
    catch
        put!(cache.cache_pool, proto)
        rethrow()
    end
end

function CommonSolve.solve!(output::RasterStack, cache::MicroMapCache)
    proto, first_I, first_result = _solve_proto_pixel!(cache)
    try
        layers = cache.problem.model.output_layers
        _write_output!(output, first_result, layers, first_I)
        return _solve_remaining!(output, cache, proto)
    catch
        put!(cache.cache_pool, proto)
        rethrow()
    end
end

# Pull worker #1 from the pool, solve the first pixel with it, return
# the worker + its result for output sizing/writing. The caller is
# responsible for returning the worker to the pool (via `_solve_remaining!`
# or its catch path).
function _solve_proto_pixel!(cache)
    first_I = first(DimIndices(cache.terrain.elevation))
    proto = take!(cache.cache_pool)
    reinit!(proto.micro, cache.init_inputs.build_inputs(proto.scratch, first_I))
    solve!(proto.micro)
    return proto, first_I, proto.micro.output
end

function _solve_remaining!(output, cache, proto)
    cache_pool = cache.cache_pool
    layers     = cache.problem.model.output_layers
    build_inputs = cache.init_inputs.build_inputs
    put!(cache_pool, proto)

    pixel_indices = DimIndices(cache.terrain.elevation)
    npixels = length(pixel_indices)
    work = Channel{eltype(pixel_indices)}(max(npixels - 1, 1))
    for (k, I) in enumerate(pixel_indices)
        k == 1 && continue
        put!(work, I)
    end
    close(work)

    nworkers = cache_pool.sz_max
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            c = take!(cache_pool)
            try
                for I in work
                    reinit!(c.micro, build_inputs(c.scratch, I))
                    solve!(c.micro)
                    _write_output!(output, c.micro.output, layers, I)
                end
            finally
                put!(cache_pool, c)
            end
        end
    end
    return output
end

"""
    CommonSolve.solve(problem::MicroMapProblem) -> RasterStack

Shortcut for `solve!(init(problem))`.
"""
CommonSolve.solve(problem::MicroMapProblem) = solve!(init(problem))

# ---------------------------------------------------------------------------
# Output stack allocation and per-pixel writes
# ---------------------------------------------------------------------------

function _allocate_output(model::MicroModel, terrain, proto, layers::Tuple, year::Integer)
    xy = dims(terrain, (X, Y))
    ti = Ti(_ti_datetime_axis(year, model.days, model.hours))
    extra = (
        soil    = (ti, Dim{:depth}(model.depths)),
        profile = (ti, Dim{:height}(model.heights)),
        scalar  = (ti,),
    )
    return RasterStack(NamedTuple(map(layers) do spec
        _layer_name(spec) => _allocate_layer(proto, spec, xy, extra)
    end))
end

# Build the Ti axis as actual DateTimes. `days` are day-of-year integers
# (mid-month for monthly mode; 1:365 for daily); `hours` are 0..23
# floats; the result is one DateTime per (day, hour) row of the inner
# Microclimate output.
function _ti_datetime_axis(year::Integer, days::AbstractVector{<:Integer},
                           hours::AbstractVector)
    out = Vector{DateTime}(undef, length(days) * length(hours))
    base = DateTime(year, 1, 1)
    k = 1
    for d in days, h in hours
        out[k] = base + Day(d - 1) + Hour(round(Int, h))
        k += 1
    end
    return out
end

function _allocate_layer(proto, spec::LayerSpec{<:Any, K}, xy, extra) where K
    T  = typeof(first(_layer_source(proto, spec)))
    ds = (xy..., extra[K]...)
    return Raster(zeros(T, map(length, ds)...), ds)
end

# `I` is a tuple of dim selectors (`(X(i), Y(j))` from `DimIndices`), so the
# spatial axes are addressed by name and the storage order of `rast` is
# irrelevant. The trailing `Ti`/extra-dim selectors target the time and
# depth/height axes we constructed in `_allocate_layer`.
@inline function _write_slice!(rast::AbstractArray{<:Any,3}, src::AbstractVector, I::Tuple)
    @inbounds for t in eachindex(src)
        rast[I..., Ti(t)] = src[t]
    end
    return nothing
end
@inline function _write_slice!(rast::AbstractArray{<:Any,4}, src::AbstractMatrix, I::Tuple)
    # All four indices wrapped as dims so DimensionalData's setindex dispatches
    # without ambiguity. The trailing dim is `Dim{:depth}` for soil layers and
    # `Dim{:height}` for profile layers — extract its base type from the raster
    # itself so this method handles both kinds.
    extra_dim = basetypeof(last(dims(rast)))
    @inbounds for d in axes(src, 2), t in axes(src, 1)
        rast[I..., Ti(t), extra_dim(d)] = src[t, d]
    end
    return nothing
end

@inline function _write_output!(output, result, layers::Tuple, I::Tuple)
    unrolled_map(layers) do spec
        _write_slice!(output[_layer_name(spec)], _layer_source(result, spec), I)
        nothing
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Show methods
# ---------------------------------------------------------------------------
#
# Each forcing has a "status" — where the value will come from at solve
# time. `_forcing_status(model)` classifies each role purely from the model;
# `_forcing_status(problem)` overlays the problem's `data` overrides on top.
# The same printer renders both.

# Per-role classifier for the bare model. Returns `(label, missing)`
# where `missing == true` means the user must supply data on the problem.
_role_status(::Val{:dem}, model) =
    _source_status(model.dem_source, "data.dem")
_role_status(::Val{:weather}, model) =
    _source_status(model.weather_source, "data.weather")
_role_status(::Val{:landcover}, model) =
    model.landcover_source === nothing ?
        (label = "(none)", missing = false) :
        (label = string(model.landcover_source), missing = false)
_role_status(::Val{:surface_albedo}, model) =
    _surface_property_status(model.surface_albedo_source, model.landcover_source,
                             :surface_albedo)
_role_status(::Val{:roughness_height}, model) =
    _surface_property_status(model.roughness_height_source, model.landcover_source,
                             :roughness_height)

_source_status(::Nothing, data_key::String) =
    (label = "<missing — pass `$data_key`>", missing = true)
_source_status(source, _) =
    (label = string(source), missing = false)

_surface_property_status(::Nothing, ::Nothing, key) =
    (label = "<missing — pass landcover_source or `data.$key`>", missing = true)
_surface_property_status(::Nothing, lc, _) =
    (label = "from $lc defaults", missing = false)
_surface_property_status(::NamedTuple, ::Nothing, key) =
    (label = "<missing landcover_source for class weights>", missing = true)
_surface_property_status(::NamedTuple, lc, _) =
    (label = "user weights via $lc", missing = false)
_surface_property_status(s::Type{<:RasterDataSources.RasterDataSource}, _, _) =
    (label = "from $s", missing = false)
_surface_property_status(::Raster, _, _) =
    (label = "user Raster", missing = false)
_surface_property_status(v, _, _) =
    (label = "constant $v", missing = false)

const _ROLES = (:dem, :weather, :landcover, :surface_albedo, :roughness_height)

# Apply problem.data overrides over the model-level role statuses. An
# override for a special key replaces the model's label with "user override".
function _role_statuses(problem::MicroMapProblem)
    data = problem.data
    map(_ROLES) do role
        status = _role_status(Val(role), problem.model)
        if haskey(data, role)
            return (role, (label = "user override (was: $(status.label))", missing = false))
        end
        return (role, status)
    end
end

# show methods

function Base.show(io::IO, ::MIME"text/plain", model::MicroMapModel)
    println(io, "MicroMapModel")
    statuses = map(role -> (role, _role_status(Val(role), model)), _ROLES)
    _print_role_table(io, statuses)
    println(io, "  lapse_rate_model:  ", model.lapse_rate_model)
    println(io, "  output_layers:     ",
            join((_layer_name(spec) for spec in model.output_layers), ", "))
    missing_roles = [string(role) for (role, s) in statuses if s.missing]
    if !isempty(missing_roles)
        println(io)
        println(io, "missing — pass via `data` on the MicroMapProblem:")
        for role in missing_roles
            println(io, "  ", role)
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", problem::MicroMapProblem)
    println(io, "MicroMapProblem")
    println(io, "  area:   ", problem.area)
    println(io, "  years:  ", problem.years)
    println(io)
    println(io, "forcings:")
    _print_role_table(io, _role_statuses(problem))
    _print_overrides(io, problem.data, "")
    init = problem.init
    if init !== nothing
        println(io)
        println(io, "init:")
        for k in keys(init)
            v = init[k]
            summary_str = v === nothing ? "nothing (solver default)" :
                          (v isa AbstractArray ? "$(eltype(v))[$(length(v))]" : repr(v))
            println(io, "  ", rpad(string(k) * ":", 20), " ", summary_str)
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", cache::MicroMapCache)
    problem = cache.problem
    nx, ny = size(cache.terrain.elevation)
    println(io, "MicroMapCache  (", nx, "×", ny, " grid, ", nx * ny, " pixels)")
    println(io, "  cache_pool:  ", cache.cache_pool.sz_max, " workers")
    println(io)
    println(io, "forcings:")
    _print_role_table(io, _role_statuses(problem))
    _print_overrides(io, problem.data, "")
end


function _print_role_table(io::IO, statuses)
    width = maximum(length(string(role)) for (role, _) in statuses)
    for (role, status) in statuses
        marker = status.missing ? "⚠ " : "  "
        println(io, "  ", marker, rpad(string(role) * ":", width + 2), " ", status.label)
    end
end

function _print_overrides(io::IO, data::NamedTuple, indent::String)
    canonical_keys = filter(k -> !_is_special_key(k), keys(data))
    if isempty(canonical_keys)
        return false
    end
    println(io, indent, "overrides (in data):")
    for k in canonical_keys
        println(io, indent, "  ", k)
    end
    return true
end
