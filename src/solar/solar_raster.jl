"""
    SolarRasterProblem

Compute solar irradiance across a terrain raster for given days and hours.

# Fields
- `dem_source`       — RasterDataSources source type (e.g. `SRTM`) or a
                       pre-loaded `Raster`.
- `extent`           — `Extent(X=(lon_min,lon_max), Y=(lat_min,lat_max))` or a
                       4-tuple `(lon_min, lon_max, lat_min, lat_max)`.
- `days`             — `AbstractVector{<:Real}` — days of year (e.g. `[172.0]`).
- `hours`            — `AbstractVector{<:Real}` — hours of day (e.g. `6.0:1.0:18.0`).
- `albedo`           — surface reflectance (default `0.15`).
- `n_horizon_angles` — number of horizon-angle directions (default `32`).
- `solar_model`      — `SolarProblem` instance (default `SolarProblem()`).
"""
@kwdef struct SolarRasterProblem
    dem_source
    extent = nothing
    days::AbstractVector{<:Real}
    hours::AbstractVector{<:Real} = 0.0:1.0:23.0
    albedo::Float64 = 0.15
    n_horizon_angles::Int = 32
    solar_model::SolarProblem = SolarProblem()
end

_normalise_solar_extent(e::Extent) = e
_normalise_solar_extent(t::NTuple{4}) =
    Extent(X = (t[1], t[2]), Y = (t[3], t[4]))

"""
    CommonSolve.solve(prob::SolarRasterProblem) -> RasterStack

Compute solar irradiance across terrain for each pixel in the DEM.

Returns a `RasterStack` with dims `(X, Y, Ti)` and layers:
- `global_terrain`     — total irradiance on the slope-facing surface (W/m²)
- `direct_horizontal`  — direct beam on horizontal surface (W/m²)
- `diffuse_horizontal` — diffuse irradiance on horizontal surface (W/m²)
- `global_horizontal`  — total irradiance on horizontal surface (W/m²)

The `Ti` lookup values are fractional day-of-year (`day + hour/24`).
Work is distributed across `Threads.nthreads()` workers — start Julia with
`julia --threads=auto` to use all available cores.
"""
function CommonSolve.solve(prob::SolarRasterProblem)
    # 1. Load DEM
    dem = if prob.dem_source isa Raster
        prob.dem_source
    else
        ext = _normalise_solar_extent(prob.extent)
        _load_dem(prob.dem_source, ext)
    end

    # 2. Compute terrain grids — slope, aspect, horizon angles, lat/lon, pressure
    terrain = compute_terrain_grids(dem; n_horizon_angles = prob.n_horizon_angles)

    days  = collect(Float64, prob.days)
    hours = collect(Float64, prob.hours)
    nmax   = prob.solar_model.wavelength_count
    ndays  = length(days)
    nsteps = ndays * length(hours)
    nx     = size(terrain.elevation, X)
    ny     = size(terrain.elevation, Y)

    # 3. Pre-extract raw parent arrays — avoids DimensionalData indexing overhead
    # in the hot pixel loop (pure array access, no lookup dispatch).
    elev_p  = parent(terrain.elevation)
    slope_p = parent(terrain.slope)
    asp_p   = parent(terrain.aspect)
    hz_p    = parent(terrain.horizon_angles)
    press_p = parent(terrain.atmospheric_pressure)
    lat_p   = parent(terrain.latitude)
    lon_p   = parent(terrain.longitude)

    # 4. Per-worker scratch pool — solar_radiation! mutates `out` and `buffers`
    # in place, so each concurrent worker needs its own copy.
    # Channel-based pool (same pattern as the microclimate solve) avoids any
    # dependence on threadid(), which can exceed nthreads() when Julia is started
    # with an interactive thread pool (e.g. --threads=16,1).
    build_scratch() = (;
        out     = allocate_output_arrays(nsteps, ndays, nmax),
        buffers = allocate_buffers(nmax, prob.solar_model.diffuse_model),
    )
    npixels  = nx * ny
    nworkers = min(Threads.nthreads(), npixels)
    proto    = build_scratch()
    scratch_pool = Channel{typeof(proto)}(nworkers)
    put!(scratch_pool, proto)
    for _ in 2:nworkers
        put!(scratch_pool, build_scratch())
    end

    # 5. Plain Float64 output — avoids Unitful broadcasting overhead during the
    # hot loop. Units are re-attached once on the final Raster construction.
    data_gt  = zeros(Float64, nx, ny, nsteps)
    data_dh  = zeros(Float64, nx, ny, nsteps)
    data_dfh = zeros(Float64, nx, ny, nsteps)
    data_gh  = zeros(Float64, nx, ny, nsteps)

    # 6. Work queue — ix column indices, workers drain concurrently.
    # Each worker then loops over all iy rows for its ix, keeping the inner
    # iy loop tight on a single scratch buffer.
    work = Channel{Int}(nx)
    for ix in 1:nx
        put!(work, ix)
    end
    close(work)

    albedo      = prob.albedo
    solar_model = prob.solar_model
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            sc = take!(scratch_pool)
            try
                for ix in work
                    for iy in 1:ny
                        st = SolarTerrain(;
                            elevation            = elev_p[ix, iy],
                            slope                = slope_p[ix, iy],
                            aspect               = asp_p[ix, iy],
                            horizon_angles       = hz_p[ix, iy],
                            albedo,
                            atmospheric_pressure = press_p[ix, iy],
                            latitude             = lat_p[ix, iy],
                            longitude            = lon_p[ix, iy],
                        )
                        solar_radiation!(sc.out, sc.buffers, solar_model;
                                         solar_terrain = st, days, hours)
                        # Scalar ustrip per timestep — no temporary allocation.
                        @inbounds for t in 1:nsteps
                            data_gt[ix, iy, t]  = ustrip(u"W/m^2", sc.out.global_terrain[t])
                            data_dh[ix, iy, t]  = ustrip(u"W/m^2", sc.out.direct_horizontal[t])
                            data_dfh[ix, iy, t] = ustrip(u"W/m^2", sc.out.diffuse_horizontal[t])
                            data_gh[ix, iy, t]  = ustrip(u"W/m^2", sc.out.global_horizontal[t])
                        end
                    end
                end
            finally
                put!(scratch_pool, sc)
            end
        end
    end

    # 7. Wrap as unitful Rasters with (X, Y, Ti) dims
    xy_dims = dims(terrain.elevation, (X, Y))
    ti_vals = vcat([d .+ hours ./ 24 for d in days]...)
    ti_dim  = Ti(ti_vals)

    return RasterStack((;
        global_terrain     = Raster(data_gt  .* u"W/m^2", (xy_dims..., ti_dim)),
        direct_horizontal  = Raster(data_dh  .* u"W/m^2", (xy_dims..., ti_dim)),
        diffuse_horizontal = Raster(data_dfh .* u"W/m^2", (xy_dims..., ti_dim)),
        global_horizontal  = Raster(data_gh  .* u"W/m^2", (xy_dims..., ti_dim)),
    ))
end
