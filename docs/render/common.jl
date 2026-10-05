# Shared setup for the render scripts. Each script runs a simulation too slow, or data too large, for the docs
# build, and writes figures and tables into docs/src/assets/generated/, which are committed.
#
#   julia --project=<docs environment> --threads=auto docs/render/<name>.jl
#
# Data are read from ENV["RASTERDATASOURCES_PATH"] (default c:/Spatial_Data/).

haskey(ENV, "RASTERDATASOURCES_PATH") || (ENV["RASTERDATASOURCES_PATH"] = "c:/Spatial_Data/")

using MicroclimateMapper, Microclimate, RasterDataSources, Rasters, Dates, Unitful, Statistics, CairoMakie
using Microclimate: example_soil_profile, example_soil_properties_model, example_soil_hydraulic_model
using Rasters.Extents: Extent
using Printf
using ZarrDatasets   # ERA5 is read from a Zarr store

include(joinpath(@__DIR__, "..", "figure_helpers.jl"))
using .FigureHelpers

const GENERATED = mkpath(joinpath(@__DIR__, "..", "src", "assets", "generated"))
const DEPTHS = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
const HEIGHTS = [0.01, 1.2]u"m"

micro_model(; snow = true, heights = HEIGHTS, kw...) = MicroModel(; depths = DEPTHS, heights,
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
    snow_model = snow ? SnowModel() : NoSnow(), kw...)

map_model(; micro = micro_model(), dem = SRTM, weather = CRUCL2, kw...) = MicroMapModel(;
    micro_model = micro, dem_source = dem, weather_source = weather,
    surface_albedo_source = 0.15, roughness_height_source = 0.004u"m", kw...)

start(; snow = true) = snow ? (; soil_moisture = fill(0.2, length(DEPTHS)), snow_depth = 0.0u"cm") :
                              (; soil_moisture = fill(0.2, length(DEPTHS)))

save_figure(name, fig) = (save(joinpath(GENERATED, name), fig); println("  wrote ", name))

function save_table(name, header, rows)
    open(joinpath(GENERATED, name), "w") do io
        println(io, join(header, ","))
        foreach(row -> println(io, join(row, ",")), rows)
    end
    println("  wrote ", name)
end

"Run `f`, returning its result and the wall time in seconds; `nothing` and the error if it fails."
function timed(label, f)
    println("── ", label)
    try
        t = @elapsed out = f()
        println("   ", round(t; digits = 1), " s")
        return out, t
    catch e
        println("   FAILED: ", sprint(showerror, e)[1:min(end, 400)])
        return nothing, NaN
    end
end

const MACHINE = string(Sys.CPU_NAME, ", ", Sys.CPU_THREADS, " threads, Julia ", VERSION, ", ", Threads.nthreads(), " used")
