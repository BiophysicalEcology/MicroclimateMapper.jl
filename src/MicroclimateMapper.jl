module MicroclimateMapper

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
using GeoInterface
using RasterDataSources
using NCDatasets       # triggers Rasters NCDatasets extension for NetCDF support
using ArchGDAL         # triggers Rasters ArchGDAL extension for GeoTIFF support
# ZarrDatasets activates Rasters' Zarr extension for ARCO-ERA5 — currently
# blocked by a JSON version conflict with the user's dev RasterDataSources
# checkout. Until that's resolved, users running ERA5 should add
# `using ZarrDatasets` themselves in their script.
using Rasters
using Rasters: X, Y, Ti, Near, Between, lookup, setcrs, crs
using Rasters.Lookups: Intervals, Center, Sampled, Regular, Irregular,
    ForwardOrdered, ReverseOrdered, span, order, sampling
using Rasters.DimensionalData: Dim, DimArray, MergedLookup, hasdim, At, dims
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
    # Weather sources
    NCEP,
    # Microclimate drivers
    MicroMapModel,
    MicroRasterProblem,
    MicroVectorProblem,
    MicroMapCache,
    solve,
    solve!,
    init,
    reinit!,
    canonical_unit,
    strip_to_canonical,
    load_cpcsoil,
    terrain,
    # Solar output
    SolarOutputLayer,
    SOLAR_BROADBAND,
    SOLAR_PAR,
    SOLAR_UVB,
    SOLAR_NIR

include("atmosphere/aerosol.jl")
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
include("climate/crucl2.jl")
include("climate/cpcsoil.jl")
include("climate/ncep.jl")
include("climate/ncep_hourly.jl")
include("climate/awap.jl")
include("climate/era5.jl")
include("climate/gridmet.jl")
include("common.jl")
include("solar/solar_output.jl")
include("raster.jl")
include("vector.jl")

end
