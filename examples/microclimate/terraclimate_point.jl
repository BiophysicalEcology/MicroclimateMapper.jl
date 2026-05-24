# Example: microclimate simulation at a point using TerraClimate monthly data
#
# This script demonstrates the two-step API:
#   1. get_weather  — download TerraClimate data and build environment objects
#   2. simulate_microclimate — assemble terrain/soil structs and run the model
#
# The location (-89.46°, 43.14°) matches the NicheMapR validation dataset so
# outputs can be compared against micro_terra.R.
#
# For a real application, terrain parameters (elevation, slope, aspect,
# horizon angles) would typically be derived from a high-resolution DEM using
# Geomorphometry.jl:
#
#   using Geomorphometry, RasterDataSources
#   dem = getraster(SRTM30, ...)
#   elevation      = Geomorphometry.extract_elevation(dem, lon, lat)
#   slope, aspect  = Geomorphometry.slope_aspect(dem)
#   horizon_angles = Geomorphometry.horizon_angles(dem, lon, lat)
#
# Below we use flat terrain at a known elevation for simplicity.

using BiophysicalGrids
using Microclimate
using SolarRadiation
using FluidProperties
using Unitful

lon, lat = -89.4557, 43.1379
elevation = 270.0u"m"

# ---------------------------------------------------------------------------
# Step 1: download and prepare monthly weather forcing
# ---------------------------------------------------------------------------
# The lapse rate correction adjusts TerraClimate grid-cell temperatures
# to the site elevation. grid_elevation defaults to 0 m (sea level) when
# not supplied — provide it for better accuracy.
weather = get_weather(TerraClimate, lon, lat;
    ystart = 2000,
    elevation,
    # grid_elevation defaults to elevation (no lapse correction).
    # Provide the WorldClim 2.5-arcmin grid elevation at this cell for lapse correction,
    # e.g. grid_elevation = 308.0u"m" if WorldClim reports that for the cell.
    vapour_pressure_method = GoffGratch(),    # swap to Teten() or Huang()
    lapse_rate_type = EnvironmentalLapseRate(), # swap to DryAdiabaticLapseRate()
)

# ---------------------------------------------------------------------------
# Step 2: construct site
# ---------------------------------------------------------------------------
# Site bundles the location, geometry, baseline pressure, surface albedo and
# roughness height — what used to be split across `SolarTerrain` and
# `MicroTerrain` in the previous Microclimate API.
site = Site(;
    latitude             = lat * u"°",
    longitude            = lon * u"°",
    elevation,
    slope                = 0.0u"°",
    aspect               = 0.0u"°",
    horizon_angles       = fill(0.0u"°", 24),
    sky_view_fraction    = 1.0,
    albedo               = 0.15,
    roughness_height     = 0.004u"m",
    atmospheric_pressure = atmospheric_pressure(elevation),
)

# ---------------------------------------------------------------------------
# Step 3: construct soil thermal and hydraulics models
# ---------------------------------------------------------------------------
# bulk_density and mineral_density now live on the hydraulics model — the
# soil properties model reads them through the energy-balance plumbing.
depths = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"

# Organic litter cap: top two nodes get organic-soil thermal properties,
# remaining nodes use the mineral defaults (NicheMapR `cap = 1` equivalent).
n = length(depths)
mineral_conductivity_vec  = fill(2.5u"W/m/K", n);  mineral_conductivity_vec[1:2]  .= 0.2u"W/m/K"
mineral_heat_capacity_vec = fill(870.0u"J/kg/K", n); mineral_heat_capacity_vec[1:2] .= 1920.0u"J/kg/K"

soil_properties_model = CampbelldeVriesSoilProperties(;
    de_vries_shape_factor   = 0.1,          # 0.33 for organic, 0.1 for mineral
    mineral_conductivity    = mineral_conductivity_vec,
    mineral_heat_capacity   = mineral_heat_capacity_vec,
    recirculation_power     = 4.0,
    return_flow_threshold   = 0.162,
)

soil_hydraulic_model = example_soil_hydraulics(depths;
    bulk_density    = 1.3u"Mg/m^3",
    mineral_density = 2.56u"Mg/m^3",
    root_density    = fill(0.0, length(depths))u"m/m^3",
)

# ---------------------------------------------------------------------------
# Step 4: optionally apply a climate change scenario
# ---------------------------------------------------------------------------
# Swap Historical for TerraClimate{Plus2C} or TerraClimate{Plus4C} to shift
# the weather forcing by the TerraClimate scenario deltas.  The baseline
# weather download (Step 1) is always from the historical record; the scenario
# only modifies the environment structs passed to simulate_microclimate.
weather_scenario = apply_climate_scenario(Historical, weather, lon, lat)
# weather_scenario = apply_climate_scenario(TerraClimate{Plus2C}, weather, lon, lat; ystart = 2000)
# weather_scenario = apply_climate_scenario(TerraClimate{Plus4C}, weather, lon, lat; ystart = 2000)

# ---------------------------------------------------------------------------
# Step 5: build the model and simulate
# ---------------------------------------------------------------------------
# Prescribed soil moisture comes from the monthly weather data, expanded to
# (ndepths × nmonths). Swap to `DynamicSoilMoisture(; ...)` to solve moisture
# at each hourly timestep instead.
precomputed_soil_moisture = let sm = get(weather_scenario, :soil_moisture_monthly, nothing)
    isnothing(sm) ? nothing : repeat(sm', length(depths), 1)
end

model = MicroModel(;
    days   = weather_scenario.days,
    hours  = collect(0.0:1.0:23.0),
    depths,
    heights = [0.01, 2.0]u"m",
    soil_properties_model,
    soil_hydraulic_model,
    vapour_pressure_equation = GoffGratch(),  # swap to Teten() or Huang()
    config = MicroConfig(;
        time_mode              = NonConsecutiveDayMode(; iterations_per_day = 10),
        convergence            = SoilTemperatureConvergenceTolerance(;
            tolerance = 0.1u"K", max_iterations_per_day = 10),
        soil_moisture_strategy = PrescribedSoilMoisture(; precomputed_soil_moisture),
    ),
)

result = simulate_microclimate(model, site, weather_scenario; clearsky = false)

# Quick inspection of outputs
soil_T = result.soil_temperature          # Matrix (timesteps × depths) in K
air_T  = result.profile.air_temperature  # Matrix (timesteps × heights) in K
@show size(soil_T)
@show result.soil_temperature[1, :]      # first hour, all depths

# # ---------------------------------------------------------------------------
# # Sensitivity: swap vapour pressure equation
# # The method must be set consistently in get_weather (humidity lapse correction),
# # apply_climate_scenario (scenario humidity adjustment), and the MicroModel
# # (boundary-layer and soil-energy calculations).
# # ---------------------------------------------------------------------------
# vp = Huang()   # faster than GoffGratch(); also try Teten() for maximum speed
# weather_huang = get_weather(TerraClimate, lon, lat;
#     ystart = 2000, elevation, vapour_pressure_method = vp,
# )
# scenario_huang = apply_climate_scenario(Historical, weather_huang, lon, lat;
#     vapour_pressure_method = vp,
# )
# model_huang = ConstructionBase.setproperties(model; vapour_pressure_equation = vp)
# result_huang = simulate_microclimate(model_huang, site, scenario_huang)

# # ---------------------------------------------------------------------------
# # Sensitivity: dry adiabatic lapse rate (only weather changes; model unchanged)
# # ---------------------------------------------------------------------------
# weather_dry = get_weather(TerraClimate, lon, lat;
#     ystart = 2000, elevation,
#     lapse_rate_type        = DryAdiabaticLapseRate(),
#     vapour_pressure_method = GoffGratch(),
# )
# scenario_dry = apply_climate_scenario(Historical, weather_dry, lon, lat;
#     vapour_pressure_method = GoffGratch(),
# )
# result_dry = simulate_microclimate(model, site, scenario_dry)

# # ---------------------------------------------------------------------------
# # Multi-year run (2000–2002) — rebuild the model with the longer day vector
# # ---------------------------------------------------------------------------
# weather_3yr = get_weather(TerraClimate, lon, lat;
#     ystart = 2000, yfinish = 2002, elevation,
#     vapour_pressure_method = GoffGratch(),
# )
# scenario_3yr = apply_climate_scenario(Historical, weather_3yr, lon, lat;
#     vapour_pressure_method = GoffGratch(),
# )
# model_3yr = ConstructionBase.setproperties(model; days = scenario_3yr.days)
# result_3yr = simulate_microclimate(model_3yr, site, scenario_3yr)
# @show size(result_3yr.soil_temperature)  # should be (36*24, ndepths)

# ---------------------------------------------------------------------------
# Visual comparisons against NicheMapR (run manually, not in CI)
# ---------------------------------------------------------------------------
# Requires Plots, CSV, DataFrames — install in your environment if needed:
#   using Pkg; Pkg.add(["Plots", "CSV", "DataFrames"])

using CSV, DataFrames, Plots

let
    data_dir = joinpath(@__DIR__, "..", "..", "test", "data", "micro_terra")
    soil   = CSV.read(joinpath(data_dir, "soil_monthly_terra.csv"),   DataFrame)
    metout = CSV.read(joinpath(data_dir, "metout_monthly_terra.csv"), DataFrame)

    n = nrow(soil)   # 288 = 12 months × 24 hours (year 2000)
    t = 1:n

    # ---- NicheMapR reference vectors ----------------------------------------
    soiltemps_nmr = Matrix(soil[:,   ["D0cm","D2.5cm","D5cm","D10cm","D15cm","D20cm","D30cm","D50cm","D100cm","D200cm"]])
    ta1cm_nmr  = (metout.TALOC .+ 273.15) .* 1u"K"
    ta2m_nmr   = (metout.TAREF .+ 273.15) .* 1u"K"
    rh1cm_nmr  =  metout.RHLOC ./ 100.0
    rh2m_nmr   =  metout.RH    ./ 100.0
    vel1cm_nmr =  metout.VLOC  .* 1u"m/s"
    vel2m_nmr  =  metout.VREF  .* 1u"m/s"

    # ---- Julia output matrices (ntimesteps × nheights) ----------------------
    air_temperature_matrix = result.profile.air_temperature
    humidity_matrix        = result.profile.relative_humidity
    wind_matrix            = result.profile.wind_speed

    depths_labels = ["$(round(ustrip(u"cm", d); digits=1)) cm"
                     for d in [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200]u"cm"]

    # ---- Soil temperatures --------------------------------------------------
    soil0_julia = ustrip.(u"°C", u"K".(result.soil_temperature[t, 1]))
    soil0_nmr   = soiltemps_nmr[:, 1]
    soil_ylim   = (floor(min(minimum(soil0_julia), minimum(soil0_nmr))) - 1,
                   ceil( max(maximum(soil0_julia), maximum(soil0_nmr))) + 1)

    p_st = plot(layout=(3, 3), size=(900, 800),
                title=reshape(depths_labels, 1, :), legend=:outertop)
    for col in 1:9
        plot!(p_st, t, ustrip.(u"°C", u"K".(result.soil_temperature[t, col]));
              sp=col, label="Julia",     color=:red,   ylabel="°C", ylims=soil_ylim)
        plot!(p_st, t, soiltemps_nmr[:, col];
              sp=col, label="NicheMapR", color=:black)
    end
    display(p_st)

    # ---- Atmospheric profiles -----------------------------------------------
    p_atm = plot(layout=(3, 2), size=(900, 800), legend=:outertop)

    plot!(p_atm, t, ustrip.(u"°C", u"K".(air_temperature_matrix[t, 1]));
          sp=1, label="Julia", color=:red, title="Air temp 1 cm", ylabel="°C")
    plot!(p_atm, t, ustrip.(u"°C", ta1cm_nmr);
          sp=1, label="NicheMapR", color=:black)

    plot!(p_atm, t, ustrip.(u"°C", u"K".(air_temperature_matrix[t, 2]));
          sp=2, label="Julia", color=:red, title="Air temp 2 m")
    plot!(p_atm, t, ustrip.(u"°C", ta2m_nmr);
          sp=2, label="NicheMapR", color=:black)

    plot!(p_atm, t, humidity_matrix[t, 1];
          sp=3, label="Julia", color=:red, title="RH 1 cm", ylabel="–")
    plot!(p_atm, t, rh1cm_nmr;
          sp=3, label="NicheMapR", color=:black)

    plot!(p_atm, t, humidity_matrix[t, 2];
          sp=4, label="Julia", color=:red, title="RH 2 m")
    plot!(p_atm, t, rh2m_nmr;
          sp=4, label="NicheMapR", color=:black)

    plot!(p_atm, t, ustrip.(u"m/s", wind_matrix[t, 1]);
          sp=5, label="Julia", color=:red, title="Wind 1 cm", ylabel="m/s")
    plot!(p_atm, t, ustrip.(u"m/s", vel1cm_nmr);
          sp=5, label="NicheMapR", color=:black)

    plot!(p_atm, t, ustrip.(u"m/s", wind_matrix[t, 2]);
          sp=6, label="Julia", color=:red, title="Wind 2 m")
    plot!(p_atm, t, ustrip.(u"m/s", vel2m_nmr);
          sp=6, label="NicheMapR", color=:black)

    display(p_atm)
end
