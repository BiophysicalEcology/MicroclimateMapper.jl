# Gridded solar radiation across a mountainous terrain using a downloaded SRTM DEM.
#
# DEM source:  SRTM via RasterDataSources (CSI-CGIAR mirror, ~90 m / 3 arcsecond)
#
# Pipeline:
#   1. Download SRTM tile and crop to ~100 × 100 pixel study area
#   2. Reproject to UTM for metric slope/aspect/horizon calculations
#   3. Compute slope, aspect, and terrain horizon angles (Geomorphometry.jl)
#   4. Run SolarRadiation.jl per pixel for each hour
#   5. Integrate daily total and plot results (CairoMakie)
#
# Extra packages beyond BiophysicalGrids core (install once):
#   using Pkg; Pkg.add(["ArchGDAL", "Plots"])

using BiophysicalGrids
using RasterDataSources
using Rasters, ArchGDAL
using SolarRadiation
using FluidProperties
using Unitful
using Printf
using Plots
import Plots: heatmap, plot, savefig

# ============================================================================
# Configuration
# ============================================================================

# Study area — change center_lon/center_lat and location_name to move to a new location.
# extent_lat/extent_lon control the N–S and E–W span (decimal degrees).
# At mid-latitudes, 0.0833° lat ≈ 0.120° lon ≈ 100 SRTM pixels (~9 km).
location_name = "Chamonix"
center_lon    = 6.87    # °E
center_lat    = 45.92   # °N
extent_lat    = 0.0833  # degrees latitude
extent_lon    = 0.120   # degrees longitude
location_tag  = lowercase(replace(location_name, " " => "_"))

region = Extent(
    X = (center_lon - extent_lon / 2, center_lon + extent_lon / 2),
    Y = (center_lat - extent_lat / 2, center_lat + extent_lat / 2),
)

simulation_day   = 172               # day of year (1–365); 172 = 21 June (summer solstice)
hour_step        = 1.0
hours_of_day     = 6.0:hour_step:20.0
default_albedo   = 0.2
n_horizon_angles = 32

# ============================================================================
# Step 1: Download SRTM tile and crop to study area
# ============================================================================

println("Downloading SRTM DEM and reprojecting to UTM...")
(; utm_dem, x_coords_utm, y_coords_utm, nx_utm, ny_utm, cs) =
    load_utm_dem(region)
println("  UTM grid: $(nx_utm) × $(ny_utm) pixels, " *
        "cell size ≈ $(round(cs[1]; digits=1)) × $(round(cs[2]; digits=1)) m")

# ============================================================================
# Steps 2: Slope, aspect, lat/lon, horizon angles, unit-tagged rasters
# ============================================================================

println("Computing terrain grids (slope, aspect, horizons)...")
terrain_grids = compute_terrain_grids(utm_dem; n_horizon_angles)
(; dem_data, data_is_xy, y_descending,
   elevation_m, slope_deg, aspect_deg,
   latitude_deg, longitude_deg, pressure_r,
   horizons_u) = terrain_grids

center_point          = Point([center_lon, center_lat])
simulation_month      = ceil(Int, (simulation_day / 365.25) * 12)  # approximate month for aerosol lookup
aerosol_optical_depth = get_aerosol_optical_depth(center_point, 0.01, simulation_month)
solar_model = SolarProblem(; scattered_uv = false, aerosol_optical_depth)

# ============================================================================
# Step 3: Compute solar radiation for each hour
# ============================================================================

println("Computing solar radiation — day $simulation_day " *
        "(hours $(first(hours_of_day))–$(last(hours_of_day)))...")

hours_vec = collect(hours_of_day)

# Returns (ny, nx, nhours) array of W/m² values (missing for no-data pixels).
# Each pixel's solar_radiation call covers all hours at once; columns parallelised.
global_terrain = solar_radiation_grid(
    terrain_grids, solar_model,
    [Float64(simulation_day)], hours_vec;
    albedo  = default_albedo,
    extract = r -> [ustrip(u"W/m^2", v) for v in r.global_terrain],
)

# ============================================================================
# Step 4: Daily-integrated radiation (trapezoidal rule)
# ============================================================================

println("Integrating daily total radiation...")
daily_Wh = zeros(ny_utm, nx_utm)

for k in 1:(length(hours_vec) - 1)
    dt = hours_vec[k + 1] - hours_vec[k]
    for j in 1:nx_utm, i in 1:ny_utm
        v1, v2 = global_terrain[i, j, k], global_terrain[i, j, k + 1]
        (!ismissing(v1) && !ismissing(v2)) && (daily_Wh[i, j] += dt * (v1 + v2) / 2)
    end
end

daily_MJ = daily_Wh .* 0.0036  # Wh/m² → MJ/m²/day
valid_d  = filter(!iszero, vec(daily_MJ))
println("  Daily range: $(round(minimum(valid_d); digits=1)) – " *
        "$(round(maximum(valid_d); digits=1)) MJ/m²/day")

# ============================================================================
# Step 5: Plot results
# ============================================================================

println("Plotting...")
fmt_h(h) = @sprintf("%02d:%02d", floor(Int, h), round(Int, (h - floor(h)) * 60))

# ascending_y(y, m) flips both y and the (ny,nx) matrix so heatmap gets ascending y.
# Slope/aspect rasters are unit-tagged; strip units and reorder to (ny, nx) first.
to_mat(r) = map(x -> ismissing(x) ? NaN : ustrip(x), r)
to_ny_nx(r) = data_is_xy ? permutedims(to_mat(r)) : to_mat(r)

y_plt, dem_plt    = ascending_y(y_coords_utm, dem_data)
_,     slope_plt  = ascending_y(y_coords_utm, to_ny_nx(slope_deg))
_,     aspect_plt = ascending_y(y_coords_utm, to_ny_nx(aspect_deg))

common_kw = (; aspect_ratio = :equal, xlabel = "Easting (m)", ylabel = "Northing (m)")

# ── Fig. 1: Terrain properties ─────────────────────────────────────────────
p_elev = heatmap(x_coords_utm, y_plt, dem_plt;
    color = :terrain, title = "Elevation (m)",
    clims = extrema(filter(!isnan, vec(dem_plt))), common_kw...)
p_slope = heatmap(x_coords_utm, y_plt, slope_plt;
    color = :YlOrRd, title = "Slope (°)",
    clims = (0.0, maximum(filter(!isnan, vec(slope_plt)))), common_kw...)
p_aspect = heatmap(x_coords_utm, y_plt, aspect_plt;
    color = :hsv, title = "Aspect (°)", clims = (0.0, 360.0), common_kw...)

display(plot(p_elev, p_slope, p_aspect;
    layout = (1, 3), size = (1200, 420), left_margin = 4Plots.mm,
    plot_title = "Terrain — $location_name  (SRTM ~90 m,  ~$(nx_utm)×$(ny_utm) pixels)"))
savefig("$(location_tag)_terrain.png")
println("  Saved $(location_tag)_terrain.png")

# ── Fig. 2: Horizon angles (4 cardinal directions) ─────────────────────────
hz_dirs  = [(1, "N 0°"), (7, "E 90°"), (13, "S 180°"), (19, "W 270°")]
horizons = ustrip.(u"°", horizons_u)  # raw Float64 for plotting
hz_valid = filter(!isnan, vec(horizons))
hz_clims = isempty(hz_valid) ? (0.0, 1.0) : extrema(hz_valid)

hz_panels = [heatmap(x_coords_utm, y_plt, ascending_y(y_coords_utm, horizons[:, :, d])[2];
    color = :YlOrRd, title = "Horizon $lbl", clims = hz_clims,
    colorbar_title = "°", common_kw...) for (d, lbl) in hz_dirs]

display(plot(hz_panels...; layout = (1, 4), size = (1400, 420), left_margin = 4Plots.mm,
    plot_title = "Terrain horizon angles — $location_name"))
savefig("$(location_tag)_horizons.png")
println("  Saved $(location_tag)_horizons.png")

# ── Fig. 3: Solar radiation — 2×2 hourly panels ────────────────────────────
all_vals  = Float64.(filter(!ismissing, vec(global_terrain)))
s_clims   = (0.0, maximum(all_vals))
panel_idx = round.(Int, range(2, length(hours_vec) - 1; length = 4))

solar_panels = [heatmap(x_coords_utm, y_plt,
    ascending_y(y_coords_utm, Float64.(coalesce.(global_terrain[:, :, pi], NaN)))[2];
    color = :inferno, title = "Hour $(fmt_h(hours_vec[pi]))",
    clims = s_clims, colorbar_title = "W/m²", common_kw...)
    for pi in panel_idx]

display(plot(solar_panels...; layout = (2, 2), size = (1100, 900), left_margin = 4Plots.mm,
    plot_title = "Solar radiation — Day $simulation_day, $location_name"))
savefig("$(location_tag)_solar_panel.png")
println("  Saved $(location_tag)_solar_panel.png")

# ── Fig. 4: Daily total radiation ──────────────────────────────────────────
display(heatmap(x_coords_utm, y_plt, ascending_y(y_coords_utm, daily_MJ)[2];
    color = :inferno, colorbar_title = "MJ/m²/day",
    title = "Daily total solar radiation — Day $simulation_day, $location_name",
    size = (850, 720), left_margin = 6Plots.mm, common_kw...))
savefig("$(location_tag)_solar_daily.png")
println("  Saved $(location_tag)_solar_daily.png")

# ── Fig. 5: Animated solar radiation — one frame per hour ──────────────────
println("Animating...")

solar_clims = (0.0, maximum(all_vals))

anim_solar = @animate for k in 1:length(hours_vec)
    frame = Float64.(coalesce.(global_terrain[:, :, k], NaN))
    _, frame_plt = ascending_y(y_coords_utm, frame)
    heatmap(x_coords_utm, y_plt, frame_plt;
        color = :inferno, clims = solar_clims,
        title = "Global terrain radiation — $(fmt_h(hours_vec[k]))\nDay $simulation_day, $location_name",
        colorbar_title = "W/m²", common_kw...,
        titlefontsize = 9, size = (700, 600),
        left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
end
gif(anim_solar, "$(location_tag)_solar.gif"; fps = 4)
println("  Saved $(location_tag)_solar.gif")

println("\nDone! $(nx_utm)×$(ny_utm) pixel grid, day $simulation_day, " *
        "$(length(hours_vec)) hours simulated.")
println("Daily solar range: $(round(minimum(valid_d); digits=1)) – " *
        "$(round(maximum(valid_d); digits=1)) MJ/m²/day")
