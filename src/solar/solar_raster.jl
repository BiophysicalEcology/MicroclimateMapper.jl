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

    # 3. Pre-allocate SolarRadiation buffers (reused across all pixels)
    days  = collect(Float64, prob.days)
    hours = collect(Float64, prob.hours)
    nmax    = prob.solar_model.wavelength_count
    ndays   = length(days)
    nsteps  = ndays * length(hours)
    out     = allocate_output_arrays(nsteps, ndays, nmax)
    buffers = allocate_buffers(nmax, prob.solar_model.diffuse_model)

    # 4. Pre-allocate output arrays as plain (X, Y, Ti) data
    W_per_m2 = typeof(0.0u"W/m^2")
    nx = size(terrain.elevation, X)
    ny = size(terrain.elevation, Y)
    data_gt  = zeros(W_per_m2, nx, ny, nsteps)
    data_dh  = zeros(W_per_m2, nx, ny, nsteps)
    data_dfh = zeros(W_per_m2, nx, ny, nsteps)
    data_gh  = zeros(W_per_m2, nx, ny, nsteps)

    # 5. Pixel loop — integer indexing avoids DimIndices assignment complexity
    albedo = prob.albedo
    for ix in 1:nx, iy in 1:ny
        st = SolarTerrain(;
            elevation            = terrain.elevation[X(ix), Y(iy)],
            slope                = terrain.slope[X(ix), Y(iy)],
            aspect               = terrain.aspect[X(ix), Y(iy)],
            horizon_angles       = terrain.horizon_angles[X(ix), Y(iy)],
            albedo,
            atmospheric_pressure = terrain.atmospheric_pressure[X(ix), Y(iy)],
            latitude             = terrain.latitude[X(ix), Y(iy)],
            longitude            = terrain.longitude[X(ix), Y(iy)],
        )
        solar_radiation!(out, buffers, prob.solar_model;
                         solar_terrain = st, days, hours)
        data_gt[ix, iy, :]  .= out.global_terrain
        data_dh[ix, iy, :]  .= out.direct_horizontal
        data_dfh[ix, iy, :] .= out.diffuse_horizontal
        data_gh[ix, iy, :]  .= out.global_horizontal
    end

    # 6. Wrap as Rasters with (X, Y, Ti) dims
    xy_dims = dims(terrain.elevation, (X, Y))
    ti_vals = vcat([d .+ hours ./ 24 for d in days]...)
    ti_dim  = Ti(ti_vals)

    global_terrain     = Raster(data_gt,  (xy_dims..., ti_dim))
    direct_horizontal  = Raster(data_dh,  (xy_dims..., ti_dim))
    diffuse_horizontal = Raster(data_dfh, (xy_dims..., ti_dim))
    global_horizontal  = Raster(data_gh,  (xy_dims..., ti_dim))

    return RasterStack((; global_terrain, direct_horizontal,
                          diffuse_horizontal, global_horizontal))
end
