# microclimate_grid.jl
#
# Take a `MicroModel`, a DEM data-source type, and a weather data-source type;
# produce a `RasterStack` of microclimate outputs over a spatial extent and
# time range.

using Rasters: Extent
using Rasters.DimensionalData: unrolled_map

# ---------------------------------------------------------------------------
# Data-source dispatch
# ---------------------------------------------------------------------------

"""
    _load_dem(::Type{<:RasterDataSources.RasterDataSource}, area::Extent) -> Raster

Load a DEM for the given lon/lat `area`. One method per supported DEM source.
"""
function _load_dem end

# Fractional buffer applied to source raster bounds so resampling onto the
# chosen template never hits an edge.
const _AREA_BUFFER = 0.05

function _load_dem(::Type{SRTM}, area::Extent)
    lon_min, lon_max = area.X
    lat_min, lat_max = area.Y
    dlon = (lon_max - lon_min) * _AREA_BUFFER
    dlat = (lat_max - lat_min) * _AREA_BUFFER
    buffered = Extent(X = (lon_min - dlon, lon_max + dlon),
                      Y = (lat_min - dlat, lat_max + dlat))
    # SRTM tiles are a 5°×5° non-overlapping grid with `missing` paths for
    # ocean tiles; cat tiles along their spatial axes — no resampling.
    tile_paths = getraster(SRTM;
        bounds = (lon_min - dlon, lat_min - dlat, lon_max + dlon, lat_max + dlat))
    any(ismissing, tile_paths) && error(
        "Some SRTM tiles are missing for area bounds " *
        "($lon_min, $lat_min, $lon_max, $lat_max); only ocean-free areas are supported."
    )
    tiles = map(tile_paths) do p
        read(crop(Raster(p; lazy = true, missingval = Int16(0));
                  to = buffered, touches = true))
    end
    # `tile_paths` indexed `[y_tile, x_tile]` per RasterDataSources convention.
    rows = [reduce((a, b) -> cat(a, b; dims = X), @view tiles[i, :])
            for i in axes(tiles, 1)]
    dem_full = reduce((a, b) -> cat(a, b; dims = Y), rows)
    return dem_full[X(Between(lon_min, lon_max)), Y(Between(lat_min, lat_max))]
end

"""
    _load_weather(::Type{<:RasterDataSources.RasterDataSource}, area, years; scenario) -> RasterStack

Load every layer needed by the microclimate model over `area` and `years`,
cropped to `area`, with a `Ti` axis spanning all months. One method per
supported weather source.
"""
function _load_weather end

_load_weather(::Type{<:TerraClimate}, area::Extent, years; scenario = Historical) =
    load_terraclimate(area, years; scenario)

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

@inline function _write_slice!(rast::AbstractArray{<:Any,3}, src::AbstractVector,
                               I::CartesianIndex{2})
    @inbounds for t in eachindex(src)
        rast[I, t] = src[t]
    end
    return nothing
end
@inline function _write_slice!(rast::AbstractArray{<:Any,4}, src::AbstractMatrix,
                               I::CartesianIndex{2})
    @inbounds for d in axes(src, 2), t in axes(src, 1)
        rast[I, t, d] = src[t, d]
    end
    return nothing
end

@inline function _write_output!(output, result, layers::Tuple, I::CartesianIndex{2})
    unrolled_map(layers) do spec
        _write_slice!(output[_layer_name(spec)], _layer_source(result, spec), I)
        nothing
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

"""
    microclimate_grid(model, dem_source, weather_source;
                      area, years, scenario = nothing,
                      init = (soil_temperature = nothing,
                              soil_moisture    = fill(0.42 * 0.25, length(model.depths))),
                      landcover_source = nothing,
                      surface_albedo   = nothing,
                      roughness_height = nothing,
                      output_layers    = _DEFAULT_OUTPUT_LAYERS,
                      lapse_rate_type  = EnvironmentalLapseRate(),
                  ) -> RasterStack

Run the microclimate model over every pixel of `area` using `dem_source`
for terrain and `weather_source` for forcing.

The model holds all simulation choices and is constant across pixels —
`init`/`reinit!` reuses one `MicroCache` across the whole grid so per-pixel
cost is just constructing fresh `MicroInputs` and a `solve!`. The output
`RasterStack` is pre-allocated once; each pixel writes into its own slice.

# Required positional arguments
- `model`: a fully-built `MicroModel`
- `dem_source`: DEM data-source type (e.g. `SRTM`)
- `weather_source`: weather data-source type (e.g. `TerraClimate{Historical}`)

# Required keyword arguments
- `area::Extents.Extent`: spatial extent with `X` (longitude) and `Y` (latitude) ranges
- `years::AbstractRange`: years of weather data to include

# Optional keyword arguments
- `scenario`: a `WarmingScenario` type (`Plus2C`, `Plus4C`) applied as a
  delta on top of the historical baseline. `nothing` (default) skips it.
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
    model::MicroModel,
    dem_source::Type{<:RasterDataSources.RasterDataSource},
    weather_source::Type{<:RasterDataSources.RasterDataSource};
    area::Extent,
    years::AbstractRange,
    scenario::Union{Nothing,Type{<:RasterDataSources.WarmingScenario}} = nothing,
    # TODO: initial conditions should come from a climatology source.
    init = (soil_temperature = nothing,
            soil_moisture    = fill(0.42 * 0.25, length(model.depths))),
    landcover_source = nothing,
    surface_albedo   = nothing,
    roughness_height = nothing,
    output_layers    = _DEFAULT_OUTPUT_LAYERS,
    lapse_rate_type::LapseRate = EnvironmentalLapseRate(),
)
    # Use the coarser-resolution source as the template and resample the
    # finer one onto it — one solve per template pixel, no duplicate work.
    dem_native = _load_dem(dem_source, area)
    weather_kw = isnothing(scenario) ? (;) : (; scenario)
    weather    = _load_weather(weather_source, area, years; weather_kw...)
    template   = first(values(weather))
    dem        = Rasters.resample(dem_native; to = template, method = :average)
    terrain    = compute_terrain_grids(dem)

    albedo_grid    = _resolve_surface_property(surface_albedo,
        landcover_source, template, area, default_landcover_albedo)
    roughness_grid = _resolve_surface_property(roughness_height,
        landcover_source, template, area, default_landcover_roughness)

    # Shared, immutable, build-once-for-the-run constants used by
    # `assemble_weather!`. Per-pixel elevation/pressure/lat/lon are swapped
    # into `flat_terrain_template` via `setproperties` inside the call.
    cloud_constants = (;
        solar_model = SolarProblem(),
        hours       = collect(0.0:1.0:23.0),
        flat_terrain_template = SolarTerrain(;
            elevation = 0.0u"m", slope = 0.0u"°", aspect = 0.0u"°",
            horizon_angles = Fill(0.0u"°", 24),
            albedo = 0.15,
            atmospheric_pressure = atmospheric_pressure(0.0u"m"),
            latitude = 0.0u"°", longitude = 0.0u"°",
        ),
    )

    ndays = length(MONTHLY_BASE_DAYS) * length(years)
    nsteps = length(cloud_constants.hours) * ndays
    nmax = cloud_constants.solar_model.wavelength_count
    allocate_scratch() = (;
        weather = allocate_weather_buffers(length(years)),
        solar = (;
            out     = allocate_output_arrays(nsteps, ndays, nmax),
            buffers = allocate_buffers(nmax, cloud_constants.solar_model.diffuse_model),
        ),
        cloud_constants,
    )

    vapour_pressure_method = model.vapour_pressure_equation

    function build_inputs(scratch, I::CartesianIndex{2})
        horizon_angles = terrain.horizon_angles[I]
        site = Site(;
            elevation            = terrain.elevation[I],
            slope                = terrain.slope[I],
            aspect               = terrain.aspect[I],
            latitude             = terrain.latitude[I],
            longitude            = terrain.longitude[I],
            atmospheric_pressure = terrain.atmospheric_pressure[I],
            horizon_angles,
            sky_view_fraction    = _sky_view_from_horizon(horizon_angles),
            albedo               = albedo_grid[I],
            roughness_height     = roughness_grid[I],
        )
        env = assemble_weather!(scratch, weather, site, I;
            vapour_pressure_method, lapse_rate_type)
        MicroInputs(; site,
            env.environment_minmax, env.environment_daily, env.environment_hourly,
            initial_soil_temperature = init.soil_temperature,
            initial_soil_moisture    = init.soil_moisture,
        )
    end

    # Solve the first pixel to size the output stack from its result eltypes.
    # That cache becomes worker #1; remaining workers get fresh ones.
    indices = vec(CartesianIndices(terrain.elevation))
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
