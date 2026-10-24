# Example: microclimate forecast at a single point from an ECMWF SEAS5
# seasonal forecast, driven via `MicroVectorProblem`.
#
# Needs a CDS API key, and the Terms of Use accepted (once) at
# https://cds.climate.copernicus.eu/datasets/seasonal-original-single-levels

using MicroclimateMapper
using Microclimate
using RasterDataSources
using Microclimate: example_microclimate_problem, example_soil_profile
using Rasters
using Unitful
using Dates
using ZarrDatasets  # for init_source = ERA5ECMWFLand
#using Plots

ENV["CDSAPI_KEY"] = "your key"
ENV["RASTERDATASOURCES_PATH"] = "your folder"

points = [(144.96, -37.81)]  # (lon, lat), Melbourne

# One forecast (start month) and ensemble member per run: member 0 is the
# control; loop `member = 0:50` for the full ensemble. A forecast covers
# 215 days from its start, so `dates` must lie within that window.
weather_source = seas5_forecast(Date(2026, 9, 1); member = 0)
dates = Date(2026, 9, 1):Day(1):Date(2026, 11, 30)

# ---------------------------------------------------------------------------
# Build the model
# ---------------------------------------------------------------------------

depths = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
           20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0] .* u"cm"
heights = [0.01, 1.2]u"m"

micro_model = MicroModel(;
    depths,
    heights,
    soil_properties_model = CampbelldeVriesSoilProperties(;
        de_vries_shape_factor = 0.1,
        recirculation_power = 4.0,
        return_flow_threshold = 0.162,
    ),
    soil_hydraulic_model  = example_soil_hydraulic_model(),
    snow_model            = SnowModel(),
    config                = MicroConfig(soil_moisture_strategy = DynamicSoilMoisture()),
)

model = MicroMapModel(;
    micro_model,
    dem_source = SRTM,
    weather_source,
    # SEAS5 has no soil state to start from: seed it from ERA5-Land at the
    # forecast start (same CDS key).
    init_source = ERA5ECMWFLand,
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
# The first run queues 7 CDS requests (one per SEAS5 variable), typically a
# few minutes; later runs in the same 1° cells reuse the downloads.
problem = MicroVectorProblem(; model, points, dates,
    soil_profile = example_soil_profile(depths),
)

@time output = solve(problem)

# plot(u"°C".(output.soil_temperature[point=1, depth=At(0.0)]))
