module BiophysicalGrids

using Dates
using Statistics: mean

using Unitful

using Microclimate
using CommonSolve: CommonSolve
using Microclimate: DEFAULT_DEPTHS
using Microclimate: example_soil_hydraulic_model, example_soil_properties_model, example_soil_profile
using SolarRadiation
using FluidProperties
using FluidProperties: GoffGratch, Teten, Huang, VapourPressureEquation

using GeoFormatTypes
using RasterDataSources
using NCDatasets       # triggers Rasters NCDatasets extension for NetCDF support
using ArchGDAL         # triggers Rasters ArchGDAL extension for GeoTIFF support
using Rasters
using Rasters: X, Y, Ti, Near, Between, lookup
using Geomorphometry
using Geomorphometry: Horn
using StaticArrays: SVector
using FillArrays: Fill
using ConstructionBase: setproperties

export
    # Lapse rate types
    LapseRate,
    EnvironmentalLapseRate,
    DryAdiabaticLapseRate,
    SaturatedAdiabaticLapseRate,
    CustomLapseRate,
    # Lapse rate functions
    lapse_rate,
    lapse_adjust_temperature,
    dewpoint_lapse_adjust,
    rh_at_temperature,
    # Cloud / radiation
    cloud_from_srad,
    # Wind
    wind_profile_adjust,
    # Data source types (re-exported from RasterDataSources)
    SRTM,
    TerraClimate,
    EarthEnv,
    LandCover,
    Historical,
    Plus2C,
    Plus4C,
    # Weather data
    get_weather,
    apply_climate_scenario,
    # Land-cover
    load_landcover,
    landcover_weighted,
    default_landcover_albedo,
    default_landcover_roughness,
    # Microclimate types (re-exported from Microclimate.jl)
    Site,
    MicroModel,
    MicroInputs,
    MicroProblem,
    MicroConfig,
    CampbelldeVriesSoilProperties,
    CampbellSoilHydraulics,
    DailyTimeseries,
    MonthlyMinMaxEnvironment,
    HourlyTimeseries,
    FixedSoilTemperatureIterations,
    SoilTemperatureConvergenceTolerance,
    NonConsecutiveDayMode,
    ConsecutiveDayMode,
    PrescribedSoilMoisture,
    DynamicSoilMoisture,
    NoSnow,
    SnowModel,
    # Microclimate helpers
    example_soil_hydraulic_model,
    example_soil_properties_model,
    example_soil_profile,
    # Simulation
    simulate_microclimate,
    microclimate_grid,
    # Terrain utilities
    get_utm_crs,
    load_utm_dem,
    compute_terrain_grids,
    compute_horizon_angles,
    # FluidProperties vapour pressure methods (not exported by FluidProperties itself)
    VapourPressureEquation,
    GoffGratch,
    Teten,
    Huang

include("terrain/terrain_utils.jl")
include("mesoclimate/mesoclimate.jl")
include("landcover/landcover.jl")
include("weatherdata/common.jl")
include("weatherdata/terraclimate.jl")
include("weatherdata/climate_scenarios.jl")
include("microclimate_grid.jl")

end
