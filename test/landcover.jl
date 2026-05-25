using Test
using BiophysicalGrids
using BiophysicalGrids: landcover_weighted
using Rasters, Rasters.Lookups
using Rasters: X, Y, RasterStack
using Unitful

# Build a small (X, Y) extent so cells are unambiguous.
xy_dims() = (X(0.0:1.0:2.0), Y(0.0:1.0:2.0))

@testset "fractional landcover_weighted" begin
    # 3×3 grid with two classes. `forest` weight 0.1, `grass` weight 0.3.
    # Fractions sum to 100 in most pixels; one row is 50/50 so weighted
    # albedo there should be the midpoint.
    forest = Raster(Float64[
        100  100   50
         50   50   25
          0    0    0
    ], xy_dims(); name = :forest)
    grass = Raster(Float64[
          0    0   50
         50   50   75
        100  100  100
    ], xy_dims(); name = :grass)
    stack = RasterStack((; forest, grass))
    weights = (forest = 0.10, grass = 0.30)

    # 3rd argument is the landcover source — ignored for fractional stacks.
    out = landcover_weighted(stack, weights, nothing)
    @test size(out) == (3, 3)
    @test out[1, 1] ≈ 0.10                  # pure forest
    @test out[3, 1] ≈ 0.30                  # pure grass
    @test out[1, 3] ≈ (50*0.10 + 50*0.30) / 100  # 50/50 mix

    # Renormalises when fractions don't sum to 100.
    @test out[2, 1] ≈ (50*0.10 + 50*0.30) / 100   # totals 100

    # Units propagate through (roughness lengths in metres).
    weights_m = (forest = 1.0u"m", grass = 0.05u"m")
    out_m = landcover_weighted(stack, weights_m, nothing)
    @test eltype(out_m) <: Quantity
    @test out_m[1, 1] === 1.0u"m"
end

@testset "categorical landcover_weighted" begin
    # Use MODIS{MCD12Q1}'s real class codes — 4 = deciduous_broadleaf_forest,
    # 12 = croplands, 17 = water_bodies — so the source's
    # `landcover_class_codes` mapping does the lookup.
    raster = Raster(Int[
         4   4  12
         4  12  12
         4  12  17
    ], xy_dims())
    weights = (
        deciduous_broadleaf_forest = 0.10,
        croplands                  = 0.20,
        water_bodies               = 0.06,
    )

    out = landcover_weighted(raster, weights, MODIS{MCD12Q1})
    @test size(out) == (3, 3)
    @test out[1, 1] ≈ 0.10
    @test out[2, 2] ≈ 0.20
    @test out[3, 3] ≈ 0.06

    # Pixels with codes outside the codes table fall back to the first weight.
    raster_unknown = Raster(Int[
         4   99  12
         4   4   12
         4   12  17
    ], xy_dims())
    out_unknown = landcover_weighted(raster_unknown, weights, MODIS{MCD12Q1})
    @test out_unknown[2, 1] ≈ 0.10  # fallback to first weight

    # Units propagate.
    weights_m = (
        deciduous_broadleaf_forest = 1.0u"m",
        croplands                  = 0.05u"m",
        water_bodies               = 0.0002u"m",
    )
    out_m = landcover_weighted(raster, weights_m, MODIS{MCD12Q1})
    @test eltype(out_m) <: Quantity
    @test out_m[3, 3] === 0.0002u"m"
end

@testset "default property tables" begin
    # Dispatch wiring: every supported land-cover source must declare both
    # an albedo and a roughness table whose keys match.
    for src in (EarthEnv{LandCover}, MODIS{MCD12Q1})
        albedo    = default_landcover_albedo(src)
        roughness = default_landcover_roughness(src)
        @test keys(albedo) == keys(roughness)
        @test all(0 .<= values(albedo) .<= 1)
        @test all(v -> v > 0u"m", values(roughness))
    end

    # Categorical sources additionally need a class→code mapping with
    # exactly the same keys.
    codes = BiophysicalGrids.landcover_class_codes(MODIS{MCD12Q1})
    @test keys(codes) == keys(default_landcover_albedo(MODIS{MCD12Q1}))
    @test all(c -> c isa Integer, values(codes))
end
