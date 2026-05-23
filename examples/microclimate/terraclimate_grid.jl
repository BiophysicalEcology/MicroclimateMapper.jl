# Gridded microclimate across the Chamonix SRTM DEM using TerraClimate forcing (July 2000).
#
# Pipeline:
#   1. Download SRTM DEM (~100×100 pixels) and reproject to UTM
#   2. Compute slope, aspect, horizon angles (Geomorphometry.jl)
#   3. Download TerraClimate weather (year 2000) at center pixel; slice to July
#   4. Per-pixel simulate_microclimate with lapse-rate-corrected weather
#   5. Plot soil surface T, air T at 1 cm, soil T at 5 cm, soil T at 20 cm
#      at 6 times of day: midnight, dawn, mid-morning, midday, mid-afternoon, dusk
#
# Each pixel runs a single July representative day (24 hours) with up to iterate_day=10
# passes (the default), stopping early when the maximum nodal temperature change is
# below convergence_tolerance=0.1 K, so near-surface soil temperatures converge.
# Rows are processed sequentially; each pixel is warm-started from the converged
# midnight soil-temperature profile of the pixel directly above it, reducing the
# iterations needed to converge.  Columns within each row are parallelised.
# The deep soil boundary condition is the annual mean temperature at each pixel
# (lapse-corrected from center elevation), computed from the full TerraClimate
# annual record.
#
# Threading: start Julia with `julia --threads auto` for best performance.
#
# Extra packages beyond BiophysicalGrids core (install once):
#   using Pkg; Pkg.add(["ArchGDAL", "Plots"])

using BiophysicalGrids
using RasterDataSources
using Rasters, ArchGDAL
using SolarRadiation
using FluidProperties
using FluidProperties: Teten, Huang, GoffGratch, VPLookupTable
using Unitful
using Statistics: median
using Printf
using Plots
import Plots: heatmap, plot, savefig

# ============================================================================
# Configuration
# ============================================================================

# Study area: above Chamonix, French Alps — same extent as grid_solar.jl
center_lon = 6.87     # °E
center_lat = 45.92    # °N
extent_lat = 0.0833   # ~100 SRTM pixels N–S
extent_lon = 0.120    # ~100 SRTM pixels E–W

lon_min = center_lon - extent_lon / 2
lon_max = center_lon + extent_lon / 2
lat_min = center_lat - extent_lat / 2
lat_max = center_lat + extent_lat / 2

year             = 2000
july             = 7          # month index
n_horizon_angles = 24

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 2.0]u"m"

# Vapour pressure formula used in humidity lapse correction and simulation.
# GoffGratch() is the most accurate but slowest; alternatives from FluidProperties:
#   Teten() — simple empirical, fastest
#   Huang()  — more accurate than Teten, faster than GoffGratch
#   VPLookupTable() — precomputed table, fastest for gridded runs
vp_method = VPLookupTable()


# Time snapshots: step index is 1-based (hour + 1)
snapshot_hours = collect(0:23)          # all 24 hours of the day
snapshot_steps = snapshot_hours .+ 1   # 1-based step index within the day
hour_labels    = [@sprintf("%02d:00", h) for h in snapshot_hours]
nhours         = length(snapshot_hours)   # 24

# Subset used for the static 2×3 panel plots (6 representative hours)
panel_ks  = [1, 7, 10, 13, 16, 19]    # indices into snapshot_hours: 00,06,09,12,15,18
panel_labels = ["Midnight", "Dawn", "Mid-morning", "Midday", "Mid-afternoon", "Dusk"]

# ============================================================================
# Steps 1–6: Load DEM, reproject to UTM, compute terrain grids
# ============================================================================

println("Downloading SRTM DEM and reprojecting to UTM...")
(; utm_dem, cs) = load_utm_dem(center_lon, center_lat, extent_lon, extent_lat)

println("Computing terrain rasters (slope, aspect, horizons)...")
(; elevation_m, slope_deg, aspect_deg,
   latitude_deg, longitude_deg, pressure_pa,
   horizon_angles_deg) = compute_terrain_grids(utm_dem; n_horizon_angles)

nx_utm = size(elevation_m, X)
ny_utm = size(elevation_m, Y)
println("  UTM grid: $(nx_utm) × $(ny_utm) pixels, " *
        "cell size ≈ $(round(cs[1]; digits=1)) × $(round(cs[2]; digits=1)) m")

# ============================================================================
# Step 7: TerraClimate weather — full year at center pixel, slice to July
#
# The full year is downloaded so that get_weather can compute the annual mean
# air temperature for use as the deep soil boundary condition. The July slice
# carries this annual mean in deep_soil_temperature (same value every month).
# ============================================================================

println("Obtaining TerraClimate weather for year $year...")

valid_elev    = filter(!ismissing, vec(parent(elevation_m)))
center_elev_u = median(ustrip.(u"m", valid_elev)) * u"m"

weather = get_weather(TerraClimate, center_lon, center_lat;
    ystart    = year,
    elevation = center_elev_u,
)

function extract_month(ws, m)
    mm = ws.environment_minmax
    ed = ws.environment_daily
    eh = ws.environment_hourly
    new_mm = MonthlyMinMaxEnvironment(;
        reference_temperature_min = mm.reference_temperature_min[[m]],
        reference_temperature_max = mm.reference_temperature_max[[m]],
        reference_wind_min        = mm.reference_wind_min[[m]],
        reference_wind_max        = mm.reference_wind_max[[m]],
        reference_humidity_min    = mm.reference_humidity_min[[m]],
        reference_humidity_max    = mm.reference_humidity_max[[m]],
        cloud_min                 = mm.cloud_min[[m]],
        cloud_max                 = mm.cloud_max[[m]],
        minima_times              = mm.minima_times,
        maxima_times              = mm.maxima_times,
    )
    new_ed = DailyTimeseries(;
        shade                 = ed.shade[[m]],
        soil_wetness          = ed.soil_wetness[[m]],
        surface_emissivity    = ed.surface_emissivity[[m]],
        cloud_emissivity      = ed.cloud_emissivity[[m]],
        rainfall              = ed.rainfall[[m]],
        deep_soil_temperature = ed.deep_soil_temperature[[m]],
        leaf_area_index       = ed.leaf_area_index[[m]],
    )
    # Slice the hourly pressure vector: each month occupies 24 consecutive entries
    h_range = ((m - 1) * 24 + 1):(m * 24)
    new_eh = HourlyTimeseries(;
        pressure              = isnothing(eh.pressure)              ? nothing : eh.pressure[h_range],
        reference_temperature = isnothing(eh.reference_temperature) ? nothing : eh.reference_temperature[h_range],
        reference_humidity    = isnothing(eh.reference_humidity)    ? nothing : eh.reference_humidity[h_range],
        reference_wind_speed  = isnothing(eh.reference_wind_speed)  ? nothing : eh.reference_wind_speed[h_range],
        global_radiation      = isnothing(eh.global_radiation)      ? nothing : eh.global_radiation[h_range],
        longwave_radiation    = isnothing(eh.longwave_radiation)     ? nothing : eh.longwave_radiation[h_range],
        cloud_cover           = isnothing(eh.cloud_cover)           ? nothing : eh.cloud_cover[h_range],
        rainfall              = isnothing(eh.rainfall)              ? nothing : eh.rainfall[h_range],
        zenith_angle          = isnothing(eh.zenith_angle)          ? nothing : eh.zenith_angle[h_range],
    )
    return merge(ws, (;
        environment_minmax    = new_mm,
        environment_daily     = new_ed,
        environment_hourly    = new_eh,
        days                  = [ws.days[m]],
        soil_moisture_monthly = ws.soil_moisture_monthly[[m]],
    ))
end

weather_july = extract_month(weather, july)

center_elev_m = round(ustrip(u"m",  center_elev_u); digits = 0)
tmin_july_C   = round(ustrip(u"°C", weather_july.environment_minmax.reference_temperature_min[1]); digits = 1)
tmax_july_C   = round(ustrip(u"°C", weather_july.environment_minmax.reference_temperature_max[1]); digits = 1)
tdeep_C       = round(ustrip(u"°C", weather_july.environment_daily.deep_soil_temperature[1]);       digits = 1)
println("  Center elevation: $center_elev_m m")
println("  July Tmin: $tmin_july_C °C,  Tmax: $tmax_july_C °C,  deep soil T: $tdeep_C °C")

# ============================================================================
# Step 8: Shared soil model and lapse correction helper
# ============================================================================

soil_thermal = CampbelldeVriesSoilProperties(;
    de_vries_shape_factor = 0.1,
    mineral_conductivity  = 2.5u"W/m/K",
    mineral_heat_capacity = 870.0u"J/kg/K",
    recirculation_power   = 4.0,
    return_flow_threshold = 0.162,
)

soil_hydraulics = example_soil_hydraulics(depths;
    bulk_density    = 1.3u"Mg/m^3",
    mineral_density = 2.56u"Mg/m^3",
    root_density    = fill(0.0, length(depths))u"m/m^3",
)

solar_model = SolarProblem()

# Lapse-correct temperature and humidity for a given elevation difference (Δz = pixel − center).
# Humidity is adjusted by conserving actual vapour pressure (dry adiabatic approximation)
# using BiophysicalGrids.rh_at_temperature.
function lapse_correct_weather(ws, elev_diff; method = vp_method)
    mm = ws.environment_minmax
    ed = ws.environment_daily
    T_min_new = lapse_adjust_temperature(mm.reference_temperature_min, elev_diff, EnvironmentalLapseRate())
    T_max_new = lapse_adjust_temperature(mm.reference_temperature_max, elev_diff, EnvironmentalLapseRate())
    new_mm = MonthlyMinMaxEnvironment(;
        reference_temperature_min = T_min_new,
        reference_temperature_max = T_max_new,
        reference_wind_min        = mm.reference_wind_min,
        reference_wind_max        = mm.reference_wind_max,
        # RH_min pairs with T_max time; RH_max pairs with T_min time
        reference_humidity_min    = rh_at_temperature(mm.reference_humidity_min, mm.reference_temperature_max, T_max_new, method),
        reference_humidity_max    = rh_at_temperature(mm.reference_humidity_max, mm.reference_temperature_min, T_min_new, method),
        cloud_min                 = mm.cloud_min,
        cloud_max                 = mm.cloud_max,
        minima_times              = mm.minima_times,
        maxima_times              = mm.maxima_times,
    )
    new_ed = DailyTimeseries(;
        shade                 = ed.shade,
        soil_wetness          = ed.soil_wetness,
        surface_emissivity    = ed.surface_emissivity,
        cloud_emissivity      = ed.cloud_emissivity,
        rainfall              = ed.rainfall,
        deep_soil_temperature = lapse_adjust_temperature(ed.deep_soil_temperature, elev_diff, EnvironmentalLapseRate()),
        leaf_area_index       = ed.leaf_area_index,
    )
    return merge(ws, (; environment_minmax = new_mm, environment_daily = new_ed))
end

# ============================================================================
# Step 9: Per-pixel microclimate simulation
# ============================================================================

println("Running per-pixel microclimate ($(nx_utm) × $(ny_utm) pixels, " *
        "$(Threads.nthreads()) thread(s))...")

# Precompute per-pixel lapse-corrected weather and initial-T values (single-threaded)
# to keep struct allocation out of the hot threaded path.
weather_per_pixel               = Array{Any}(undef, nx_utm, ny_utm)
mean_air_temperature_per_pixel  = fill(NaN * u"K", nx_utm, ny_utm)
deep_soil_temperature_per_pixel = fill(NaN * u"K", nx_utm, ny_utm)

for j in 1:ny_utm, i in 1:nx_utm
    elevation = elevation_m[X = i, Y = j]
    ismissing(elevation) && continue
    weather_pixel = lapse_correct_weather(weather_july, elevation - center_elev_u)
    weather_per_pixel[i, j]              = weather_pixel
    mean_air_temperature_per_pixel[i, j] = (weather_pixel.environment_minmax.reference_temperature_min[1] +
                                            weather_pixel.environment_minmax.reference_temperature_max[1]) / 2
    deep_soil_temperature_per_pixel[i, j] = weather_pixel.environment_daily.deep_soil_temperature[1]
end

# Per-pixel hourly output buffers — plain arrays during the threaded loop,
# wrapped as Rasters once the simulation finishes.
soil_temperature_surface_data = fill(NaN, nx_utm, ny_utm, nhours)
air_temperature_1cm_data      = fill(NaN, nx_utm, ny_utm, nhours)
air_temperature_2m_data       = fill(NaN, nx_utm, ny_utm, nhours)
soil_temperature_10cm_data    = fill(NaN, nx_utm, ny_utm, nhours)

# Spatial warm-start: store the converged midnight soil-temperature profile so
# adjacent rows can reuse it as initial_soil_temperature. Rows sequential,
# columns parallel.
converged_midnight_profile = fill!(Matrix{Any}(undef, nx_utm, ny_utm), nothing)
pixels_done                = Threads.Atomic{Int}(0)
n_total                    = nx_utm * ny_utm

@time for j in 1:ny_utm
    Threads.@threads :static for i in 1:nx_utm
        elevation = elevation_m[X = i, Y = j]
        latitude  = latitude_deg[X = i, Y = j]
        longitude = longitude_deg[X = i, Y = j]
        slope     = slope_deg[X = i, Y = j]
        aspect    = aspect_deg[X = i, Y = j]
        pressure  = pressure_pa[X = i, Y = j]

        if ismissing(elevation) || ismissing(latitude) || ismissing(longitude) ||
           ismissing(slope)     || ismissing(aspect)   || ismissing(pressure)
            continue
        end

        weather_pixel             = weather_per_pixel[i, j]
        mean_air_temperature_july = mean_air_temperature_per_pixel[i, j]
        deep_soil_temperature     = deep_soil_temperature_per_pixel[i, j]

        initial_soil_temperature = if j > 1 && !isnothing(converged_midnight_profile[i, j - 1])
            converged_midnight_profile[i, j - 1]
        else
            T_initial      = Vector{typeof(1.0u"K")}(undef, 10)
            T_initial[1:8] .= mean_air_temperature_july
            T_initial[9]   = (mean_air_temperature_july + deep_soil_temperature) / 2
            T_initial[10]  = deep_soil_temperature
            T_initial
        end

        site = Site(;
            latitude,
            longitude,
            elevation,
            slope,
            aspect,
            horizon_angles       = parent(horizon_angles_deg[X = i, Y = j]),
            sky_view_fraction    = 1.0,
            albedo               = 0.15,
            roughness_height     = 0.004u"m",
            atmospheric_pressure = pressure,
        )

        result = simulate_microclimate(
            site, soil_thermal, soil_hydraulics, weather_pixel;
            depths, heights, solar_model,
            initial_soil_temperature,
            vapour_pressure_equation = vp_method,
            iterate_day              = 5,
        )

        converged_midnight_profile[i, j] = collect(result.soil_temperature[1, :])

        air_temperature = result.profile.air_temperature
        @inbounds for (k, step) in enumerate(snapshot_steps)
            soil_temperature_surface_data[i, j, k] = ustrip(u"°C", result.soil_temperature[step, 1])
            air_temperature_1cm_data[i,      j, k] = ustrip(u"°C", air_temperature[step, 1])
            air_temperature_2m_data[i,       j, k] = ustrip(u"°C", air_temperature[step, 2])
            soil_temperature_10cm_data[i,    j, k] = ustrip(u"°C", result.soil_temperature[step, 4])
        end

        done = Threads.atomic_add!(pixels_done, 1)
        if Threads.threadid() == 1 && (done % nx_utm == 0 || done == n_total - 1)
            percent = round(100 * (done + 1) / n_total; digits = 1)
            print("  row $j/$ny_utm  ($(done+1)/$n_total, $percent%)   \r")
        end
    end
end
println("\nSimulation complete.")

# Wrap per-pixel outputs as 3-D Rasters with (X, Y, Ti) dimensions.
time_dim    = Ti(snapshot_hours .* u"hr")
output_dims = (dims(elevation_m, X), dims(elevation_m, Y), time_dim)
output_crs  = crs(elevation_m)

soil_temperature_surface = Raster(soil_temperature_surface_data, output_dims; crs = output_crs)
air_temperature_1cm      = Raster(air_temperature_1cm_data,      output_dims; crs = output_crs)
air_temperature_2m       = Raster(air_temperature_2m_data,       output_dims; crs = output_crs)
soil_temperature_10cm    = Raster(soil_temperature_10cm_data,    output_dims; crs = output_crs)

# ============================================================================
# Step 10: Plots — one 2×3 figure per variable
# ============================================================================

function plot_variable(raster, variable_label, filename)
    valid_values = filter(!isnan, vec(parent(raster)))
    isempty(valid_values) && return
    color_limits = (minimum(valid_values), maximum(valid_values))
    n_frames     = size(raster, Ti)
    indices      = n_frames <= 6 ? (1:n_frames) : panel_ks
    labels       = n_frames <= 6 ? hour_labels[1:n_frames] : panel_labels

    panels = [Plots.heatmap(view(raster, Ti = indices[n]);
        color = cgrad(:RdYlBu, rev = true), clims = color_limits,
        title = labels[n], colorbar_title = "°C",
        titlefontsize = 9, aspect_ratio = :equal,
        xlabel = "Easting (m)", ylabel = "Northing (m)") for n in eachindex(indices)]

    display(Plots.plot(panels...; layout = (2, 3), size = (1400, 900),
        left_margin = 5Plots.mm,
        plot_title = "$variable_label — Chamonix, July $year"))
    Plots.savefig(filename)
    println("  Saved $filename")
end

function animate_variable(raster, variable_label, filename; framerate = 4)
    valid_values = filter(!isnan, vec(parent(raster)))
    isempty(valid_values) && return
    color_limits = (minimum(valid_values), maximum(valid_values))
    n_frames     = size(raster, Ti)
    frame_labels = n_frames == length(snapshot_hours) ?
        [@sprintf("%02d:00", snapshot_hours[k]) for k in 1:n_frames] :
        [@sprintf("frame %d", k) for k in 1:n_frames]

    animation = @animate for k in 1:n_frames
        Plots.heatmap(view(raster, Ti = k);
            color = cgrad(:RdYlBu, rev = true), clims = color_limits,
            title = "$variable_label\n$(frame_labels[k]) — Chamonix, July $year",
            xlabel = "Easting (m)", ylabel = "Northing (m)",
            colorbar_title = "°C", aspect_ratio = :equal,
            titlefontsize = 9, size = (700, 600),
            left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
    end
    Plots.gif(animation, filename; fps = framerate)
    println("  Saved $filename")
end

println("Plotting...")
plot_variable(air_temperature_2m,       "Air temperature at 2 m (°C)",    "chamonix_july_air_temperature_2m.png")
plot_variable(air_temperature_1cm,      "Air temperature at 1 cm (°C)",   "chamonix_july_air_temperature_1cm.png")
plot_variable(soil_temperature_surface, "Soil surface temperature (°C)",  "chamonix_july_soil_temperature_surface.png")
plot_variable(soil_temperature_10cm,    "Soil temperature at 10 cm (°C)", "chamonix_july_soil_temperature_10cm.png")

println("Animating...")
animate_variable(air_temperature_2m,       "Air temperature at 2 m (°C)",    "chamonix_july_air_temperature_2m.gif")
animate_variable(air_temperature_1cm,      "Air temperature at 1 cm (°C)",   "chamonix_july_air_temperature_1cm.gif")
animate_variable(soil_temperature_surface, "Soil surface temperature (°C)",  "chamonix_july_soil_temperature_surface.gif")
animate_variable(soil_temperature_10cm,    "Soil temperature at 10 cm (°C)", "chamonix_july_soil_temperature_10cm.gif")

println("\nDone. $(nx_utm)×$(ny_utm) pixel grid, July $year, " *
        "$(nhours) time snapshots ($(join(snapshot_hours, ", ")) h).")
