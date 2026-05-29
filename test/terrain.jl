using Test
using MicroclimateMapper
using MicroclimateMapper: compute_terrain_grids
using Rasters, Rasters.Lookups
using Rasters: X, Y
using GeoFormatTypes
import Proj   # activates Rasters' ProjExt so cellarea / reproject work
using Unitful

# Manhattan / UTM-ish grid: 100 m cells, "Projected" CRS WKT for a generic
# UTM zone. We only need slope/aspect/horizon math, not a particular zone.
const _TEST_WKT = """PROJCS["test_utm",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1]]"""

function _synthetic_dem(values::AbstractMatrix)
    nx, ny = size(values)
    xs = X(range(500_000.0; step = 100.0, length = nx);
        sampling = Intervals(Start()), crs = GeoFormatTypes.WellKnownText(GeoFormatTypes.CRS(), _TEST_WKT))
    ys = Y(range(4_700_000.0; step = 100.0, length = ny);
        sampling = Intervals(Start()), crs = GeoFormatTypes.WellKnownText(GeoFormatTypes.CRS(), _TEST_WKT))
    return Raster(values, (xs, ys); crs = GeoFormatTypes.WellKnownText(GeoFormatTypes.CRS(), _TEST_WKT))
end

@testset "compute_terrain_grids: flat DEM" begin
    dem = _synthetic_dem(fill(Int16(100), 5, 5))
    grids = compute_terrain_grids(dem; n_horizon_angles = 8)

    @test propertynames(grids) ==
        (:elevation, :slope, :aspect, :latitude, :longitude,
         :atmospheric_pressure, :horizon_angles)
    @test size(grids.elevation) == (5, 5)
    @test unit(eltype(grids.elevation)) == u"m"
    @test unit(eltype(grids.slope)) == u"°"
    @test unit(eltype(grids.aspect)) == u"°"
    @test unit(eltype(grids.latitude)) == u"°"
    @test unit(eltype(grids.longitude)) == u"°"
    @test unit(eltype(grids.atmospheric_pressure)) == u"Pa"

    # Flat terrain → zero slope everywhere away from the edges (edge cells
    # may pick up boundary effects depending on Geomorphometry's `Horn`).
    interior_slopes = grids.slope[X(2:4), Y(2:4)]
    @test all(s -> isapprox(ustrip(u"°", s), 0; atol = 1e-3), interior_slopes)

    # Horizon angle vectors have the requested length.
    @test length(grids.horizon_angles[1, 1]) == 8
    @test eltype(grids.horizon_angles[1, 1]) <: Quantity
end

@testset "compute_terrain_grids: sloped DEM" begin
    # Plane increasing 10 m per cell to the east (X) → ~5.7° slope facing east.
    nx, ny = 6, 6
    values = Int16[10 * (i - 1) for i in 1:nx, _ in 1:ny]
    dem = _synthetic_dem(values)
    grids = compute_terrain_grids(dem; n_horizon_angles = 4)

    interior_slopes = ustrip.(u"°", grids.slope[X(2:nx-1), Y(2:ny-1)])
    expected_deg = atand(10 / 100)
    @test all(s -> isapprox(s, expected_deg; atol = 0.5), interior_slopes)
end
