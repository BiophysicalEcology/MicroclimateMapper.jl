# Example: microclimate simulation at a single point using ERA5 hourly
# reanalysis data, driven via `MicroVectorProblem`.
#
# The location (-89.46°, 43.14°) matches the NicheMapR validation dataset so
# outputs can be compared against `micro_era5` output (see the commented
# block at the bottom).
#
# Extra packages for the validation plot (install once):
#   using Pkg; Pkg.add(["Plots", "CSV", "DataFrames"])

using MicroclimateMapper
using Microclimate: example_microclimate_problem
using RasterDataSources
using Rasters
using Unitful
using Dates
using Zarr  # reads ERA5's ARCO Zarr stores
using GeoInterface: Wrappers as GIW

point = GIW.Point((-89.4557, 43.1379))
dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)

# ---------------------------------------------------------------------------
# Choose an ERA5 source
# ---------------------------------------------------------------------------
# - `ERA5{ECMWF,ARCO}`: ECMWF's Zarr store, chunked for long series at a
#   point: fast, no queue.
# - `ERA5{ECMWF,ARCO,Land}`: ERA5-Land at 0.1°, with cloud cover from ERA5.
# - `ERA5{ECMWF,OnDemand}` and `ERA5{ECMWF,OnDemand,Land}`: monthly NetCDF
#   downloads from the Climate Data Store, queued, so slower.
# - `ERA5`: Google's public store, no account, but chunked by hour over the
#   whole globe, so slow for a long series at a point.
#
# ECMWF sources need a CDS token:
# ENV["RASTERDATASOURCES_ERA5_CDS_TOKEN"] = "your token"
weather_source = ERA5{ECMWF,ARCO}

# ---------------------------------------------------------------------------
# Build the model
# ---------------------------------------------------------------------------
# The inner micro_model comes from Microclimate.jl's example builder — it
# carries sensible defaults for soil profile, hydraulics, and config.
inner = example_microclimate_problem().model

model = MicroMapModel(;
    micro_model = inner,
    dem_source = SRTM,
    weather_source,
    surface_albedo_source = 0.15,
    roughness_height_source = 0.004u"m",
)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
# ERA5 has no soil moisture here, so start from a uniform profile.
init = (; soil_moisture = fill(0.2, length(inner.soil.depths)))
problem = MicroVectorProblem(; model, points = [point], dates, init)
output = solve(problem)

@show size(output.soil_temperature)
@show output.soil_temperature[1, 1, :]   # first point, first hour, all depths

# ---------------------------------------------------------------------------
# Aerosol optical depth (optional)
# ---------------------------------------------------------------------------
# `get_aerosol_optical_depth` returns the spectral AOD vector at the
# requested location, RH fraction, and month. Wiring it into the inner
# `SolarProblem` requires rebuilding `inner.config.solar_model` — left as
# an exercise.
aod = get_aerosol_optical_depth(point, 0.01, 6)
@show length(aod)

# ---------------------------------------------------------------------------
# Visual comparison against NicheMapR (uncomment to run; requires
# Plots, CSV, DataFrames)
# ---------------------------------------------------------------------------
# using CSV, DataFrames, Plots
#
# let
#     data_dir = joinpath(@__DIR__, "..", "..", "test", "data", "micro_era5")
#     soil = CSV.read(joinpath(data_dir, "soil_era5.csv"), DataFrame)
#     metout = CSV.read(joinpath(data_dir, "metout_era5.csv"), DataFrame)
#
#     n = nrow(soil)
#     t = 1:n
#
#     # ---- NicheMapR reference vectors ----------------------------------------
#     soiltemps_nmr = Matrix(soil[:, ["D0cm","D2.5cm","D5cm","D10cm","D15cm","D20cm","D30cm","D50cm","D100cm","D200cm"]])
#     ta1cm_nmr = (metout.TALOC .+ 273.15) .* 1u"K"
#     ta2m_nmr = (metout.TAREF .+ 273.15) .* 1u"K"
#
#     # ---- Julia output (first point) -----------------------------------------
#     soil_T = view(output.soil_temperature, 1, :, :)   # (Ti, depth)
#     air_T = view(output.air_temperature, 1, :, :)   # (Ti, height)
#     depths_labels = ["$(round(ustrip(u"cm", d); digits=1)) cm"
#                      for d in [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200]u"cm"]
#
#     p_st = plot(layout=(3, 3), size=(900, 800),
#                 title=reshape(depths_labels, 1, :), legend=:outertop)
#     for col in 1:9
#         plot!(p_st, t, ustrip.(u"°C", u"K".(soil_T[t, col]));
#               sp=col, label="Julia", color=:red, ylabel="°C")
#         plot!(p_st, t, soiltemps_nmr[:, col];
#               sp=col, label="NicheMapR", color=:black)
#     end
#     display(p_st)
#
#     p_air = plot(layout=(1, 2), size=(900, 400), legend=:outertop)
#     plot!(p_air, t, ustrip.(u"°C", u"K".(air_T[t, 1]));
#           sp=1, label="Julia", color=:red, title="Air temp 1 cm", ylabel="°C")
#     plot!(p_air, t, ustrip.(u"°C", ta1cm_nmr);
#           sp=1, label="NicheMapR", color=:black)
#     plot!(p_air, t, ustrip.(u"°C", u"K".(air_T[t, 2]));
#           sp=2, label="Julia", color=:red, title="Air temp 2 m")
#     plot!(p_air, t, ustrip.(u"°C", ta2m_nmr);
#           sp=2, label="NicheMapR", color=:black)
#     display(p_air)
# end
