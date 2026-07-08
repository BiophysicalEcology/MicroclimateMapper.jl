using MicroclimateMapper
using Microclimate
using Microclimate: example_soil_hydraulic_model
using Rasters, RasterDataSources
using Rasters.Extents: Extent
using Dates, Unitful

ENV["RASTERDATASOURCES_PATH"] = "Z:"

depths = Microclimate.DEFAULT_DEPTHS
site = geocode("Alice Springs, Australia")
lon, lat = site.lon, site.lat
area = Extent(X = (lon - 0.05, lon + 0.05), Y = (lat - 0.05, lat + 0.05))

point_result = build_soil_profile(SoilGrids, (lon, lat); depths)
sp_soilgrids_point = point_result.soil_profile
campbell_b = point_result.campbell_b
air_entry_potential = point_result.air_entry_potential
saturated_conductivity = point_result.saturated_conductivity
field_capacity = point_result.field_capacity
wilting_point = point_result.wilting_point

sp_soilgrids_area = build_soil_profile(SoilGrids, area; depths).soil_profile
sp_slga = build_soil_profile(SLGA, area; depths).soil_profile

for model in (CosbyUnivariate(), CosbyMultivariate(), Campbell1985())
    sp = build_soil_profile(SoilGrids, area; depths, pedotransfer_model = model).soil_profile
    @info "$(nameof(typeof(model))): b=$(round.(sp.hydraulics.campbell_b_parameter; digits = 2))"
end

for (label, sp) in ((:soilgrids_point, sp_soilgrids_point), (:soilgrids_area, sp_soilgrids_area), (:slga, sp_slga))
    h = sp.hydraulics
    @assert all(2 .<= h.campbell_b_parameter .<= 16) "$label: campbell_b out of range"
    @assert all(0u"J/kg" .< h.air_entry_water_potential .< 20u"J/kg") "$label: air_entry_potential out of range"
    @assert all(1e-8u"kg*s/m^3" .< h.saturated_hydraulic_conductivity .< 1u"kg*s/m^3") "$label: Ksat out of range"
    @assert all(0.0u"Mg/m^3" .< sp.bulk_density .< 2.0u"Mg/m^3") "$label: bulk_density out of range"
end
@info "campbell_b" round.(campbell_b; digits = 2)
@info "air_entry_potential" round.(air_entry_potential; digits = 3)
@info "saturated_conductivity" saturated_conductivity
@info "field_capacity" round.(field_capacity; digits = 3)
@info "wilting_point" round.(wilting_point; digits = 3)

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths,
        heights = [0.01, 1.2]u"m",
        soil_properties_model = Microclimate.example_soil_properties_model(),
        soil_hydraulic_model = example_soil_hydraulic_model(),
        snow_model = NoSnow(),
    ),
    dem_source = CRUCL2,
    weather_source = CRUCL2,
    soil_moisture_source = CPCSoil,
    surface_albedo_source = 0.15,
    roughness_height_source = 0.004u"m",
    compute_terrain = false,
)

problem = MicroVectorProblem(;
    model,
    points = [site],
    dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31),
    soil_profile = sp_soilgrids_point,
    init = (; soil_moisture = fill(0.2, length(depths))),
)

@time output = solve(problem);
@info "solved ok" size(output.soil_temperature)
