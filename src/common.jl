# micro_common.jl
#
# Shared types and scaffolding for the two declarative microclimate-simulation
# entry points:
#
#   MicroRasterProblem (grid mode, src/micro_raster.jl)
#   MicroVectorProblem (points mode, src/micro_vector.jl)
#
# Both produce a `MicroMapCache` and run through the same `solve!` /
# `_solve_proto_pixel!` / `_solve_remaining!` loop — the only difference is
# the spatial layout of the forcings (a 2-D `(X, Y)` grid vs a 1-D
# `Dim{:point}` lookup carrying point coordinates).
#
# The mode-specific `init` methods load and spatially-extract their forcings
# (weather, terrain, albedo, roughness, canonical overrides) into matching
# `Raster` shapes, then hand off to the shared `_build_inputs_and_pool`
# factory which produces the per-pixel `build_inputs` closure and worker pool
# used by `solve!`.

# ---------------------------------------------------------------------------
# Output layer specs and per-pixel write
# ---------------------------------------------------------------------------

# A `LayerSpec{Name, Kind}` encodes both the result-field name and its storage
# kind in its type, so each unrolled iteration specialises on a concrete type.
# `Kind` is one of:
#   :soil    — result.<Name> is a (Ti, Depth) matrix
#   :profile — result.profile.<Name> is a (Ti, Height) matrix
#   :scalar  — result.<Name> is a (Ti,) vector
# Number of azimuth directions used for horizon-angle computation. A
# power of two makes the cardinal directions land on sample points.
const N_HORIZON_ANGLES = 32

struct LayerSpec{Name, Kind} end
LayerSpec(name::Symbol, kind::Symbol) = LayerSpec{name, kind}()

@inline _layer_name(::LayerSpec{N}) where N = N
@inline _layer_source(result, ::LayerSpec{N, :profile}) where N = getproperty(result.profile, N)
@inline _layer_source(result, ::LayerSpec{N}) where N = getproperty(result, N)

const _DEFAULT_OUTPUT_LAYERS = (
    LayerSpec(:soil_temperature, :soil),
    LayerSpec(:soil_moisture, :soil),
    LayerSpec(:air_temperature, :profile),
    LayerSpec(:relative_humidity, :profile),
    LayerSpec(:wind_speed, :profile),
    LayerSpec(:surface_water, :scalar),
    LayerSpec(:global_radiation, :scalar),
    LayerSpec(:sky_temperature, :scalar),
    LayerSpec(:snow_depth, :scalar),
)

"""
    canonical_unit(name::Symbol)

Canonical `Unitful` unit each default output layer is reported in for
storage and downstream analysis. Used by `strip_to_canonical` when
preparing a stack for unit-free I/O (e.g. NetCDF).
"""
canonical_unit(name::Symbol)              = canonical_unit(Val(name))
canonical_unit(::Val{:soil_temperature})  = u"°C"
canonical_unit(::Val{:air_temperature})   = u"°C"
canonical_unit(::Val{:sky_temperature})   = u"°C"
canonical_unit(::Val{:snow_depth})        = u"cm"
canonical_unit(::Val{:wind_speed})        = u"m/s"
canonical_unit(::Val{:surface_water})     = u"kg/m^2"
canonical_unit(::Val{:global_radiation})  = u"W/m^2"
canonical_unit(::Val{:soil_moisture})     = u"m^3/m^3"
canonical_unit(::Val{:relative_humidity}) = u"percent"

"""
    strip_to_canonical(stack::RasterStack)

Convert every layer of `stack` to its `canonical_unit` and strip the
Unitful unit, returning a Float-valued `RasterStack` suitable for writing
to a format that does not carry Unitful metadata (e.g. NetCDF).
"""
function strip_to_canonical(stack)
    names = propertynames(stack)
    return RasterStack(NamedTuple{names}(map(names) do n
        ustrip.(canonical_unit(n), stack[n])
    end))
end

# ---------------------------------------------------------------------------
# Shared types
# ---------------------------------------------------------------------------

"""
    MicroMapModel(; micro_model, dem_source, weather_source,
                    landcover_source=nothing,
                    surface_albedo_source=nothing,
                    roughness_height_source=nothing,
                    output_layers=_DEFAULT_OUTPUT_LAYERS,
                    lapse_rate_model=EnvironmentalLapseRate())

Declarative description of a microclimate run. Constant across spatial
extents and time ranges — pair with a `MicroRasterProblem` (grid) or
`MicroVectorProblem` (points) to actually run.

- `micro_model::MicroModel` — inner per-pixel physics model
- `dem_source` — DEM data-source type (e.g. `SRTM`)
- `weather_source` — weather data-source type (e.g. `TerraClimate{Historical}`)
- `landcover_source` — optional land-cover data-source type
  (e.g. `EarthEnv{LandCover}`)
- `surface_albedo_source`, `roughness_height_source` — accepted forms:
    * `nothing` (use `landcover_source`'s default class→value table)
    * a class→value `NamedTuple` (overrides the dataset default)
    * a scalar (broadcast)
    * a `Raster` (resampled/extracted to the run template)
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
    landcover_source::LCS = nothing
    surface_albedo_source::SAS = nothing
    roughness_height_source::RHS = nothing
    output_layers::OL = _DEFAULT_OUTPUT_LAYERS
    lapse_rate_model::LRT = EnvironmentalLapseRate()
end

# Per-run workspace built by `init(problem)`. Holds the spatially-extracted
# forcings, pre-resolved surface grids, and a Channel pool of per-worker
# inner `MicroCache`s + scratch buffers. The output `RasterStack` is
# allocated fresh on every `solve!` call — it does not live on the cache.
#
# Used by both `MicroRasterProblem` and `MicroVectorProblem`: only the shape of
# the spatial dim differs (`(X, Y)` vs `(Dim{:point},)`); every solve-time
# code path indexes via `I::Tuple` of dim wrappers from `DimIndices`, so the
# loop body is mode-agnostic.
mutable struct MicroMapCache{P,W,T,A,R,CO,M,POOL,SC,CC}
    problem::P                       # MicroRasterProblem/MicroVectorProblem
    weather::W                       # RasterStack
    terrain::T                       # RasterStack
    albedo_grid::A                   # Raster
    roughness_grid::R                # Raster
    canonical_overrides::CO          # NamedTuple of pre-extracted Rasters
    mask::M                          # nothing, or a 2-D Bool Raster (X, Y)
    cache_pool::POOL                 # Channel{(; micro::MicroCache, scratch)}
    init_inputs::SC                  # Resolved initial conditions + build_inputs closure
    cloud_constants::CC              # Shared SolarRadiation constants
end

# ---------------------------------------------------------------------------
# Init helpers — shared by both modes
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

# Resolve initial conditions. Every key defaults to `nothing`, meaning
# "fall back to the weather source / data override / inner-solver default"
# — never fabricated. User-supplied keys override; missing keys are filled
# from `_DEFAULT_INIT` so the build_inputs closure always sees a complete
# NamedTuple.
const _DEFAULT_INIT = (
    soil_temperature = nothing,
    soil_moisture = nothing,
    snow_depth = nothing,
    snow_temperature = nothing,
    snow_density = nothing,
)
_resolve_init(it::NamedTuple, _) = merge(_DEFAULT_INIT, it)
_resolve_init(::Nothing, _) = _DEFAULT_INIT

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
        "or supply `data.soil_moisture` on the problem."
    )
    return fill(buffer[1], length(depths))
end

# Source-loading helpers (override beats loader). Mode-agnostic: both grid
# and points init build their own `area::Extent` (rectangular for grid;
# bounding box of points for points mode) and call these.
_resolve_dem(data::NamedTuple, source, area) =
    haskey(data, :dem) ? data.dem : _load_dem(source, area)

_resolve_weather(data::NamedTuple, source, area, years) =
    haskey(data, :weather) ? data.weather : _load_weather(source, area, years)

# Normalise `area` to an Extent for loaders; geometry is kept separately for the mask.
_to_extent(area::Extent) = area
_to_extent(area) = GeoInterface.extent(area)

# `Extent` area = no mask (run every pixel); geometry = rasterise into a Bool mask.
_build_area_mask(::Extent, _template) = nothing
_build_area_mask(geom, template) = Rasters.boolmask(geom; to = template)

@inline _is_active(_I::Tuple, ::Nothing) = true
@inline _is_active(I::Tuple, mask) = mask[I...]

_first_active_index(rast, ::Nothing) = first(DimIndices(rast))
function _first_active_index(rast, mask)
    for I in DimIndices(rast)
        mask[I...] && return I
    end
    error("Area mask excludes every pixel of the template — nothing to solve.")
end

# ---------------------------------------------------------------------------
# Native surface-property resolution
# ---------------------------------------------------------------------------
#
# `_resolve_surface_native` returns a *native-resolution* representation of a
# surface property (albedo / roughness) — either a scalar (for the "broadcast
# everywhere" case) or a 2-D `(X, Y)` `Raster`. Mode-specific code then runs
# this through the appropriate spatial extractor (`_resample(...; to=template)`
# in grid mode, per-point `Near`-lookup in points mode).
#
# The forms accepted by `MicroMapModel.surface_albedo_source` /
# `.roughness_height_source` are identical to the legacy single-step path:
#   - `nothing` → dispatch landcover-source default class→value table
#   - `NamedTuple` of weights → load landcover, weight at native res
#   - scalar → return the scalar; mode wraps it in a `Fill`-Raster
#   - `Raster` → return as-is
#   - `Type{<:RasterDataSource}` → load via `load_surface_property`

# Scalar / unrecognised value: return as-is. Mode extractor builds the Fill.
_resolve_surface_native(value, _landcover_source, _area, _default_fn) = value
# User-supplied `Raster`: return as-is. Mode extractor resamples/extracts.
_resolve_surface_native(r::Raster, _landcover_source, _area, _default_fn) = r
# `nothing` → dispatch the default class→value NamedTuple for the
# landcover source, then re-enter via the NamedTuple branch.
_resolve_surface_native(::Nothing, landcover_source, area, default_fn) =
    _resolve_surface_native(default_fn(landcover_source), landcover_source, area, default_fn)
_resolve_surface_native(::Nothing, ::Nothing, _area, default_fn) =
    error("`$(default_fn)`: no value supplied. Pass a `landcover_source` " *
          "(e.g. `EarthEnv{LandCover}`, `MODIS{MCD12Q1}`), a `Type{<:RasterDataSource}` " *
          "for the property itself, or a scalar / Raster / NamedTuple.")
# Class→value NamedTuple: load landcover via `landcover_source` and weight
# at native resolution. Returned as a `Raster` carrying `commondims(landcover, (X, Y))`.
function _resolve_surface_native(
    weights::NamedTuple,
    landcover_source::Type{<:RasterDataSources.RasterDataSource},
    area, _default_fn,
)
    landcover = load_landcover(landcover_source, area)
    weighted_native = landcover_weighted(landcover, weights, landcover_source)
    spatial_dims = commondims(landcover, (X(), Y()))
    return _wrap_native_raster(weighted_native, spatial_dims)
end
# `Type{<:RasterDataSource}` as the property *source* itself (not landcover):
# load a 2-D Raster from that dataset; ignores `landcover_source`. Each
# dataset that supports this implements `load_surface_property(source, area)`.
_resolve_surface_native(
    property_source::Type{<:RasterDataSources.RasterDataSource},
    _landcover_source, area, _default_fn,
) = load_surface_property(property_source, area)

@inline _wrap_native_raster(values::Raster, _spatial_dims) = values
@inline _wrap_native_raster(values, spatial_dims) = Raster(values, spatial_dims)

# ---------------------------------------------------------------------------
# Build the per-pixel inputs closure + worker pool
# ---------------------------------------------------------------------------

# Build the cloud-derivation constants. Immutable, build-once-per-run.
function _build_cloud_constants()
    return (;
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
end

# Build the `build_inputs(scratch, I)` closure and the per-worker cache pool
# shared by both grid and points init. All spatial forcings must already be
# laid out so that `forcing[I...]` works for every `I` from
# `DimIndices(terrain.elevation)` (i.e. `(X(i), Y(j))` in grid mode and
# `(Dim{:point}(p),)` in points mode).
function _build_inputs_and_pool(;
    model, weather_source, weather, terrain,
    albedo_grid, roughness_grid, canonical_overrides,
    init_inputs, soil_moisture_available, years, cloud_constants,
)
    (; micro_model, lapse_rate_model) = model
    vapour_pressure_method = micro_model.vapour_pressure_equation

    resolution = temporal_resolution(weather_source)
    ndays = steps_per_year(resolution) * length(years)
    days = _days_of_year(resolution, length(years))
    time_mode = _time_mode(resolution)
    nsteps = length(cloud_constants.hours) * ndays
    nmax = cloud_constants.solar_model.wavelength_count
    allocate_scratch() = (;
        weather = allocate_weather_buffers(weather_source, length(years)),
        solar = (;
            out = allocate_output_arrays(nsteps, ndays, nmax),
            buffers = allocate_buffers(nmax, cloud_constants.solar_model.diffuse_model),
        ),
        cloud_constants,
    )

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
        # Snow init: fall back to MicroInputs' own defaults when the user
        # didn't supply a value. `initial_snow_density = nothing` is the
        # MicroInputs default ("use the snow model's snow_density") — pass
        # through unchanged.
        MicroInputs(; site,
            env.environment_minmax, env.environment_daily, env.environment_hourly,
            initial_soil_temperature = init_inputs.soil_temperature,
            initial_soil_moisture,
            initial_snow_depth = something(init_inputs.snow_depth, 0.0u"cm"),
            initial_snow_temperature = something(init_inputs.snow_temperature, u"K"(0.0u"°C")),
            initial_snow_density = init_inputs.snow_density,
        )
    end

    npixels = length(terrain.elevation)
    first_I = first(DimIndices(terrain.elevation))
    build_cache() = let scratch = allocate_scratch()
        (micro = CommonSolve.init(MicroProblem(micro_model, build_inputs(scratch, first_I); days, time_mode)),
         scratch)
    end

    nworkers = min(Threads.nthreads(), npixels)
    proto = build_cache()
    cache_pool = Channel{typeof(proto)}(nworkers)
    put!(cache_pool, proto)
    for _ in 2:nworkers
        put!(cache_pool, build_cache())
    end

    return (build_inputs, cache_pool)
end

# ---------------------------------------------------------------------------
# CommonSolve.solve! / solve
# ---------------------------------------------------------------------------

"""
    CommonSolve.solve!(cache::MicroMapCache) -> RasterStack
    CommonSolve.solve!(output::RasterStack, cache::MicroMapCache) -> RasterStack

Run the per-pixel microclimate solve across every spatial index of
`cache.terrain.elevation`. In grid mode that's every `(X, Y)` pixel; in
points mode it's every `Dim{:point}` entry. Work is dispatched over the
cache pool with one `Threads.@spawn` worker per pool slot.

The one-arg form allocates a fresh `RasterStack` to write into — sized
from the first pixel's result so the output's element types are concrete.
The two-arg form writes into the supplied `output` and returns it, for
callers that want to reuse storage across runs.
"""
function CommonSolve.solve!(cache::MicroMapCache)
    proto, first_I, first_result = _solve_proto_pixel!(cache)
    try
        layers = cache.problem.model.output_layers
        output = _allocate_output(cache.problem.model.micro_model, proto.micro.problem.days,
            cache.terrain, first_result, layers, first(cache.problem.years), cache.mask)
        _write_output!(output, first_result, layers, first_I)
        return _solve_remaining!(output, cache, proto, first_I)
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
        return _solve_remaining!(output, cache, proto, first_I)
    catch
        put!(cache.cache_pool, proto)
        rethrow()
    end
end

# Pull worker #1 from the pool, solve the first active pixel with it, return
# the worker + its result for output sizing/writing. The caller is
# responsible for returning the worker to the pool (via `_solve_remaining!`
# or its catch path).
function _solve_proto_pixel!(cache)
    first_I = _first_active_index(cache.terrain.elevation, cache.mask)
    proto = take!(cache.cache_pool)
    reinit!(proto.micro, cache.init_inputs.build_inputs(proto.scratch, first_I))
    solve!(proto.micro)
    return proto, first_I, proto.micro.output
end

function _solve_remaining!(output, cache, proto, first_I)
    cache_pool = cache.cache_pool
    layers = cache.problem.model.output_layers
    build_inputs = cache.init_inputs.build_inputs
    mask = cache.mask
    put!(cache_pool, proto)

    pixel_indices = DimIndices(cache.terrain.elevation)
    npixels = length(pixel_indices)
    work = Channel{eltype(pixel_indices)}(max(npixels - 1, 1))
    for I in pixel_indices
        I == first_I && continue
        _is_active(I, mask) || continue
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

# ---------------------------------------------------------------------------
# Output stack allocation and per-pixel writes
# ---------------------------------------------------------------------------

function _allocate_output(model::MicroModel, days, terrain, proto, layers::Tuple, year::Integer, mask)
    spatial_dims = dims(terrain.elevation)
    ti = Ti(_ti_datetime_axis(year, days, model.hours))
    # Dim lookups must be plain numbers, not Unitful — Unitful in dims
    # breaks NetCDF I/O and surprises selectors. Strip both to metres so
    # `depth` and `height` share a single linear axis.
    extra = (
        soil = (ti, Dim{:depth}(ustrip.(u"m", model.depths))),
        profile = (ti, Dim{:height}(ustrip.(u"m", model.heights))),
        scalar = (ti,),
    )
    return RasterStack(NamedTuple(map(layers) do spec
        _layer_name(spec) => _allocate_layer(proto, spec, spatial_dims, extra, mask)
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

function _allocate_layer(proto, spec::LayerSpec{<:Any, K}, spatial_dims, extra, mask) where K
    T = typeof(first(_layer_source(proto, spec)))
    ds = (spatial_dims..., extra[K]...)
    return _allocate_layer_array(T, ds, mask)
end

_allocate_layer_array(T, ds, ::Nothing) = Raster(zeros(T, map(length, ds)...), ds)
_allocate_layer_array(T, ds, _mask) =
    Rasters.create(nothing, T, ds; missingval = missing, fill = missing)

# `I` is a tuple of dim selectors from `DimIndices`: `(X(i), Y(j))` in grid
# mode, `(Dim{:point}(p),)` in points mode. The spatial axes are addressed
# by name so the storage order of `rast` is irrelevant. The trailing
# `Ti`/extra-dim selectors target the time and depth/height axes we
# constructed in `_allocate_layer`. Dispatch is on the **source's** rank
# (Vector = scalar, Matrix = soil/profile) so the same methods serve grid
# (rast 3-D / 4-D) and points (rast 2-D / 3-D) mode.
@inline function _write_slice!(rast, src::AbstractVector, I::Tuple)
    @inbounds for t in eachindex(src)
        rast[I..., Ti(t)] = src[t]
    end
    return nothing
end
@inline function _write_slice!(rast, src::AbstractMatrix, I::Tuple)
    # The trailing dim is `Dim{:depth}` for soil layers and `Dim{:height}`
    # for profile layers — extract its base type from the raster itself so
    # this method handles both kinds.
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
# Show helpers (shared)
# ---------------------------------------------------------------------------

# Each forcing has an origin describing where its value comes from at solve
# time. `_forcing_origin` classifies the origin from the model alone;
# `_forcing_origins(problem)` overlays the problem's `data` overrides on top.
_forcing_origin(::Val{:dem}, model) = _source_origin(model.dem_source, "data.dem")
_forcing_origin(::Val{:weather}, model) = _source_origin(model.weather_source, "data.weather")
_forcing_origin(::Val{:landcover}, model) =
    model.landcover_source === nothing ?
        (label = "(none)", missing = false) :
        (label = string(model.landcover_source), missing = false)
_forcing_origin(::Val{:surface_albedo}, model) =
    _surface_property_origin(model.surface_albedo_source, model.landcover_source, :surface_albedo)
_forcing_origin(::Val{:roughness_height}, model) =
    _surface_property_origin(model.roughness_height_source, model.landcover_source, :roughness_height)

# Apply problem.data overrides over the model-level origins. An override
# for a special key replaces the model's label with "user override".
function _forcing_origins(problem)
    data = problem.data
    map(_FORCINGS) do forcing
        origin = _forcing_origin(Val(forcing), problem.model)
        if haskey(data, forcing)
            return (forcing, (label = "user override (was: $(origin.label))", missing = false))
        end
        return (forcing, origin)
    end
end

_source_origin(::Nothing, data_key::String) =
    (label = "<missing — pass `$data_key`>", missing = true)
_source_origin(source, _) = (label = string(source), missing = false)

_surface_property_origin(::Nothing, ::Nothing, key) =
    (label = "<missing — pass landcover_source or `data.$key`>", missing = true)
_surface_property_origin(::Nothing, lc, _) =
    (label = "from $lc defaults", missing = false)
_surface_property_origin(::NamedTuple, ::Nothing, key) =
    (label = "<missing landcover_source for class weights>", missing = true)
_surface_property_origin(::NamedTuple, lc, _) =
    (label = "user weights via $lc", missing = false)
_surface_property_origin(s::Type{<:RasterDataSources.RasterDataSource}, _, _) =
    (label = "from $s", missing = false)
_surface_property_origin(::Raster, _, _) =
    (label = "user Raster", missing = false)
_surface_property_origin(v, _, _) =
    (label = "constant $v", missing = false)

const _FORCINGS = (:dem, :weather, :landcover, :surface_albedo, :roughness_height)

function Base.show(io::IO, ::MIME"text/plain", model::MicroMapModel)
    println(io, "MicroMapModel")
    origins = map(forcing -> (forcing, _forcing_origin(Val(forcing), model)), _FORCINGS)
    _print_forcing_table(io, origins)
    println(io, "  lapse_rate_model:  ", model.lapse_rate_model)
    println(io, "  output_layers:     ",
            join((_layer_name(spec) for spec in model.output_layers), ", "))
    missing_forcings = [string(forcing) for (forcing, o) in origins if o.missing]
    if !isempty(missing_forcings)
        println(io)
        println(io, "missing — pass via `data` on the problem:")
        for forcing in missing_forcings
            println(io, "  ", forcing)
        end
    end
end

function _print_forcing_table(io::IO, origins)
    width = maximum(length(string(forcing)) for (forcing, _) in origins)
    for (forcing, origin) in origins
        marker = origin.missing ? "⚠ " : "  "
        println(io, "  ", marker, rpad(string(forcing) * ":", width + 2), " ", origin.label)
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

function _print_problem_init(io::IO, problem)
    init = problem.init
    init === nothing && return
    println(io)
    println(io, "init:")
    for k in keys(init)
        v = init[k]
        summary_str = v === nothing ? "nothing (solver default)" :
                      (v isa AbstractArray ? "$(eltype(v))[$(length(v))]" : repr(v))
        println(io, "  ", rpad(string(k) * ":", 20), " ", summary_str)
    end
end
