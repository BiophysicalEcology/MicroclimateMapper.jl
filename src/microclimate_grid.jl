# microclimate_grid.jl
#
# Take a `MicroModel`, a DEM data-source type, and a weather data-source type;
# produce a `RasterStack` of microclimate outputs over a spatial extent and
# time range.

# `_load_dem` is declared in `terrain/terrain_utils.jl`; SRTM (and any other
# DEM source) contributes methods from `terrain/<source>.jl`.
#
# `_load_weather` and `assemble_weather!` are defined in climate/weather.jl;
# each source's bindings (e.g. climate/terraclimate.jl, climate/chelsa.jl,
# climate/worldclim.jl) just declare `weather_loader`, `weather_variables`,
# and optionally `primary_layers` / `fallback_layers` / `fallback_source` /
# `weather_derivations`.

# ---------------------------------------------------------------------------
# Output layer specs and per-pixel write
# ---------------------------------------------------------------------------

# A `LayerSpec{Name, Kind}` encodes both the result-field name and its storage
# kind in its type, so each unrolled iteration specialises on a concrete type.
# `Kind` is one of:
#   :soil    — `(X, Y, Ti, Depth)`,  sourced from `result.<Name>`
#   :profile — `(X, Y, Ti, Height)`, sourced from `result.profile.<Name>`
#   :scalar  — `(X, Y, Ti)`,         sourced from `result.<Name>`
struct LayerSpec{Name, Kind} end
LayerSpec(name::Symbol, kind::Symbol) = LayerSpec{name, kind}()

@inline _layer_name(::LayerSpec{N})       where N = N
@inline _layer_source(result, ::LayerSpec{N, :profile}) where N = getproperty(result.profile, N)
@inline _layer_source(result, ::LayerSpec{N})           where N = getproperty(result, N)

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
# Orchestrator
# ---------------------------------------------------------------------------

"""
    microclimate_grid(model;
                      dem_source, weather_source, area, years,
                      landcover_source = nothing,
                      init = (soil_temperature = nothing,
                              soil_moisture    = fill(0.42 * 0.25, length(model.depths))),
                      surface_albedo   = nothing,
                      roughness_height = nothing,
                      output_layers    = _DEFAULT_OUTPUT_LAYERS,
                      lapse_rate_type  = EnvironmentalLapseRate(),
                  ) -> RasterStack

Run the microclimate model over every pixel of `area` using `dem_source`
for terrain and `weather_source` for forcing. A warming scenario is
selected via the source type parameter, e.g. `TerraClimate{Plus2C}`.

The model holds all simulation choices and is constant across pixels —
`init`/`reinit!` reuses one `MicroCache` across the whole grid so per-pixel
cost is just constructing fresh `MicroInputs` and a `solve!`. The output
`RasterStack` is pre-allocated once; each pixel writes into its own slice.

# Positional arguments
- `model`: a fully-built `MicroModel`

# Required keyword arguments
- `dem_source`: DEM data-source type (e.g. `SRTM`)
- `weather_source`: weather data-source type (e.g. `TerraClimate{Historical}`)
- `area::Extents.Extent`: spatial extent with `X` (longitude) and `Y` (latitude) ranges
- `years::AbstractRange`: years of weather data to include

# Optional keyword arguments
- `init`: NamedTuple of initial conditions (`soil_temperature`,
  `soil_moisture`). TODO: pull these from a climatology source.
- `landcover_source`: a `RasterDataSources` land-cover dataset type
  (e.g. `EarthEnv{LandCover}`), or `nothing`.
- `surface_albedo`, `roughness_height`: per-pixel surface properties.
  May be `nothing` (use the dataset's `default_landcover_albedo` /
  `default_landcover_roughness` — requires `landcover_source`), a class→
  value `NamedTuple` (overrides the dataset default), a scalar
  (broadcast), or a `Raster` (resampled to the weather template).
- `output_layers`: tuple of `LayerSpec`s controlling which output rasters
  are materialised.
- `lapse_rate_type`: atmospheric lapse-rate model for elevation-correcting
  weather data from the coarser grid down to DEM-resolution pixels.
"""
function microclimate_grid(
    model::MicroModel;
    dem_source::Type{<:RasterDataSources.RasterDataSource},
    weather_source::Type{<:RasterDataSources.RasterDataSource},
    landcover_source = nothing,
    area::Extent,
    years::AbstractRange,
    # TODO: initial conditions should come from a climatology source.
    init = (soil_temperature = nothing,
            soil_moisture = fill(0.42 * 0.25, length(model.depths))),
    # TODO: how to better group these
    surface_albedo = nothing,
    roughness_height = nothing,
    output_layers = _DEFAULT_OUTPUT_LAYERS,
    # TODO: this is a model component, should we have a composite model?
    lapse_rate_type::LapseRate = EnvironmentalLapseRate(),
)
    # Use the coarser-resolution source as the template and resample the
    # finer one onto it — one solve per template pixel, no duplicate work.
    dem_native = _load_dem(dem_source, area)
    weather = _load_weather(weather_source, area, years)
    template = first(values(weather))
    dem = Rasters.resample(dem_native; to = template, method = :average)
    terrain = compute_terrain_grids(dem)

    albedo_grid = _resolve_surface_property(surface_albedo,
        landcover_source, template, area, default_landcover_albedo)
    roughness_grid = _resolve_surface_property(roughness_height,
        landcover_source, template, area, default_landcover_roughness)

    # Shared, immutable, build-once-for-the-run constants used by
    # `assemble_weather!`. Per-pixel elevation/pressure/lat/lon are swapped
    # into `flat_terrain_template` via `setproperties` inside the call.
    cloud_constants = (;
        solar_model = SolarProblem(),
        hours = collect(0.0:1.0:23.0),
        flat_terrain_template = SolarTerrain(;
            elevation = 0.0u"m", slope = 0.0u"°", aspect = 0.0u"°",
            horizon_angles = Fill(0.0u"°", 24),
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
    nmax = cloud_constants.solar_model.wavelength_count
    allocate_scratch() = (;
        weather = allocate_weather_buffers(weather_source, length(years)),
        solar = (;
            out = allocate_output_arrays(nsteps, ndays, nmax),
            buffers = allocate_buffers(nmax, cloud_constants.solar_model.diffuse_model),
        ),
        cloud_constants,
    )

    vapour_pressure_method = model.vapour_pressure_equation

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
            vapour_pressure_method, lapse_rate_type)
        MicroInputs(; site,
            env.environment_minmax, env.environment_daily, env.environment_hourly,
            initial_soil_temperature = init.soil_temperature,
            initial_soil_moisture = init.soil_moisture,
        )
    end

    # Solve the first pixel to size the output stack from its result eltypes.
    # That cache becomes worker #1; remaining workers get fresh ones.
    # `DimIndices` gives a tuple `(X(i), Y(j))` per pixel so every downstream
    # access can stay dimension-order-agnostic via dim wrappers.
    indices = vec(DimIndices(terrain.elevation))
    first_I = first(indices)
    build_cache() = let scratch = allocate_scratch()
        (micro = CommonSolve.init(MicroProblem(model, build_inputs(scratch, first_I))),
         scratch)
    end
    proto  = build_cache()
    solve!(proto.micro)
    output = _allocate_output(model, terrain, proto.micro.output, output_layers)
    _write_output!(output, proto.micro.output, output_layers, first_I)

    # Per-worker caches in a Channel pool. `init` allocates fresh buffers,
    # state and integrator per worker — the model (immutable) is shared by
    # reference, so no deepcopy is needed.
    nworkers = min(Threads.nthreads(), length(indices) - 1)
    cache_pool = Channel{typeof(proto)}(nworkers)
    put!(cache_pool, proto)
    for _ in 2:nworkers
        put!(cache_pool, build_cache())
    end
    work = Channel{eltype(indices)}(length(indices) - 1)
    for I in @view indices[2:end]
        put!(work, I)
    end
    close(work)
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            c = take!(cache_pool)
            try
                for I in work
                    reinit!(c.micro, build_inputs(c.scratch, I))
                    solve!(c.micro)
                    _write_output!(output, c.micro.output, output_layers, I)
                end
            finally
                put!(cache_pool, c)
            end
        end
    end
    return output
end

function _allocate_output(model::MicroModel, terrain, proto, layers::Tuple)
    xy = dims(terrain, (X, Y))
    ti = Ti(1:length(model.days) * length(model.hours))
    extra = (
        soil    = (ti, Dim{:depth}(model.depths)),
        profile = (ti, Dim{:height}(model.heights)),
        scalar  = (ti,),
    )
    return RasterStack(NamedTuple(map(layers) do spec
        _layer_name(spec) => _allocate_layer(proto, spec, xy, extra)
    end))
end

function _allocate_layer(proto, spec::LayerSpec{<:Any, K}, xy, extra) where K
    T = typeof(first(_layer_source(proto, spec)))
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
