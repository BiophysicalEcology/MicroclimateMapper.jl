# Terrain: Saba. The LiDAR DEM of Geomorphometry.jl and SolarRadiation.jl's tutorial, placed at 17.63 N, 63.23 W.
include(joinpath(@__DIR__, "common.jl"))
using Downloads
using Rasters.Lookups: Sampled, Intervals, Center

const YEAR = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)
const LAT0, LON0 = 17.63, -63.23
saba_file = joinpath(mkpath(joinpath(ENV["RASTERDATASOURCES_PATH"], "saba")), "saba.tif")
isfile(saba_file) || Downloads.download("https://github.com/Deltares/Geomorphometry.jl/releases/download/v0.6.0/saba.tif", saba_file)

# The DEM is in metres on a local grid with no location. Give it longitude and latitude about the island's
# centre, every `step`th cell, sea set to missing.
function saba_lonlat(step; window = nothing)
    dtm = Raster(saba_file)
    xs, ys = collect(lookup(dtm, X)), collect(lookup(dtm, Y))
    ix, iy = window === nothing ? (1:step:length(xs), 1:step:length(ys)) : (window[1][1]:step:window[1][2], window[2][1]:step:window[2][2])
    xc, yc = (first(xs) + last(xs)) / 2, (first(ys) + last(ys)) / 2
    lon = LON0 .+ (xs[ix] .- xc) ./ (111_320 * cosd(LAT0))
    lat = LAT0 .+ (ys[iy] .- yc) ./ 110_574
    vals = [ismissing(v) ? missing : Float64(v) for v in parent(dtm)[ix, iy]]
    # Ranges give regular lookups, which the package needs for its grid arithmetic.
    lon = range(first(lon), last(lon); length = length(lon))
    lat = range(first(lat), last(lat); length = length(lat))
    Raster(vals, (X(Sampled(lon; sampling = Intervals(Center()))), Y(Sampled(lat; sampling = Intervals(Center()))));
        crs = EPSG(4326), missingval = missing)
end

saba_problem(dem; kw...) = MicroRasterProblem(;
    model = map_model(; micro = micro_model(; snow = false), weather = TerraClimate{Historical}, kw...),
    area = Rasters.Extents.extent(dem), template = dem, dates = YEAR,
    soil_profile = example_soil_profile(DEPTHS), init = start(; snow = false), data = (; dem = coalesce.(dem, 0.0)))

rows = Any[]
island = saba_lonlat(8)
println("island at 40 m: ", size(island), ", land cells ", count(!ismissing, island))
out, t = timed("Saba island, 40 m", () -> solve(saba_problem(island)))
if out !== nothing
    push!(rows, ("island", "40 m", count(!ismissing, island), round(t; digits = 1)))
    times = collect(lookup(out, Ti))
    noon(m) = findfirst(t -> month(t) == m && hour(t) == 12, times)
    land(r) = Raster(ifelse.(ismissing.(parent(island)), NaN, parent(r)), dims(island))
    save_figure("saba_island.png", map_figure((
        (land(island), "Elevation (m)", :terrain, "m"),
        (land(celsius.(out.soil_temperature[depth = 1, Ti = noon(12)])), "Soil surface, December noon (°C)", :thermal, "°C"),
        (land(celsius.(out.soil_temperature[depth = 1, Ti = noon(6)])), "Soil surface, June noon (°C)", :thermal, "°C"),
        (land(celsius.(mean(out.soil_temperature[depth = 4]; dims = Ti)[Ti = 1])), "Mean soil at 10 cm (°C)", :thermal, "°C"),
    ); size = (880, 760)))
end

# Mount Scenery, the summit, at 10 m.
dtm = Raster(saba_file)
peak = Tuple(argmax(coalesce.(parent(dtm), -Inf)))
win = ((peak[1] - 60, peak[1] + 60), (peak[2] - 60, peak[2] + 60))
summit = saba_lonlat(2; window = win)
out, t = timed("Mount Scenery, 10 m", () -> solve(saba_problem(summit)))
if out !== nothing
    push!(rows, ("Mount Scenery", "10 m", count(!ismissing, summit), round(t; digits = 1)))
    times = collect(lookup(out, Ti))
    idx = [findfirst(t -> month(t) == m && hour(t) == h, times) for (m, h) in ((12, 9), (12, 15), (6, 9), (6, 15))]
    save_figure("saba_summit_drape.png", drape_panels(celsius.(out.soil_temperature[depth = 1]),
        rebuild(summit; data = Float64.(coalesce.(parent(summit), 0.0))), idx,
        ["December 09:00", "December 15:00", "June 09:00", "June 15:00"];
        title = "Soil surface temperature, Mount Scenery", label = "°C", ncols = 2))
end

# Photosynthetically active radiation alone, every 20 m.
fine = saba_lonlat(4)
par, t = timed("Saba PAR, 20 m, solar only", () -> solve(saba_problem(fine; solar_only = true, solar_output_layers = (SOLAR_PAR,))))
if par !== nothing
    push!(rows, ("island, solar only", "20 m", count(!ismissing, fine), round(t; digits = 1)))
    daily_par = sum(par.par; dims = Ti)[Ti = 1] .* u"hr"
    land20(r) = Raster(ifelse.(ismissing.(parent(fine)), NaN, ustrip.(u"MJ/m^2", parent(r)) ./ 12), dims(fine))
    save_figure("saba_par.png", map_figure(((land20(daily_par), "PAR, mean of the 12 days (MJ m⁻² d⁻¹)", :viridis, "MJ m⁻² d⁻¹"),); ncols = 1, size = (520, 440)))
end
save_table("saba.csv", ["area", "resolution", "land cells", "seconds"], rows)
println("machine: ", MACHINE)
