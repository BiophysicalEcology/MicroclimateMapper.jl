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
# ZarrDatasets activates Rasters' Zarr extension for ARCO-ERA5 — currently
# blocked by a JSON version conflict with the user's dev RasterDataSources
# checkout. Until that's resolved, users running ERA5 should add
# `using ZarrDatasets` themselves in their script.
using Rasters
using Rasters: X, Y, Ti, Near, Between, lookup
using Rasters.DimensionalData: unrolled_map, basetypeof
using Rasters.Extents
using Rasters.Extents: Extent
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
    # Data source types (re-exported from RasterDataSources)
    SRTM,
    TerraClimate,
    CHELSA,
    WorldClim,
    NCEP,
    SurfaceGauss,
    AWAP,
    ERA5,
    Climate,
    Future,
    CMIP5,
    CMIP6,
    RCP26, RCP45, RCP60, RCP85,
    SSP126, SSP245, SSP370, SSP585,
    EarthEnv,
    LandCover,
    MODIS,
    MCD12Q1,
    Historical,
    Plus2C,
    Plus4C,
    # Weather data
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
    DailyMinMaxEnvironment,
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
    microclimate_grid

include("terrain/terrain_utils.jl")
include("terrain/srtm.jl")
include("mesoclimate/lapse_rate.jl")
include("mesoclimate/cloud.jl")
include("landcover/landcover.jl")
include("landcover/earthenv.jl")
include("landcover/modis.jl")
include("climate/weather.jl")
include("climate/terraclimate.jl")
include("climate/chelsa.jl")
include("climate/worldclim.jl")
include("climate/ncep.jl")
include("climate/awap.jl")
include("climate/era5.jl")
include("microclimate_grid.jl")

end
