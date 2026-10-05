# Maps: Palm Springs and Mount San Jacinto at two resolutions, and Mont Aigoual draped over its terrain.
include(joinpath(@__DIR__, "common.jl"))

const YEAR = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)   # CRUCL2: 12 representative days
const PALM_SPRINGS = (-116.545, 33.830)
const SAN_JACINTO = (-116.679, 33.814)
const PS_SJ = Extent(X = (-116.85, -116.40), Y = (33.65, 34.00))
const AIGOUAL = Extent(X = (3.554, 3.608), Y = (44.095, 44.149))

srtm(area, factor) = (r = read(crop(Raster(SRTM; extent = area, lazy = true); to = area, touches = true));
                      factor == 1 ? r : aggregate(mean, r, factor))

grid(area, template) = MicroRasterProblem(; model = map_model(), area, template,
    dates = YEAR, soil_profile = example_soil_profile(DEPTHS), init = start())

hottest(out) = celsius.(maximum(out.soil_temperature[depth = 1]; dims = Ti)[Ti = 1])
coldest(out) = celsius.(minimum(out.soil_temperature[depth = 1]; dims = Ti)[Ti = 1])
deepest_snow(out) = ustrip.(u"cm", maximum(out.snow_depth; dims = Ti)[Ti = 1])

rows = Any[]
sites = ("Palm Springs" => PALM_SPRINGS, "San Jacinto" => SAN_JACINTO)
results = Dict()
if isempty(ARGS) || "pssj" in ARGS
for (label, factor) in (("1.8 km", 20), ("900 m", 10))
    template = srtm(PS_SJ, factor)
    out, t = timed("Palm Springs–San Jacinto at $label", () -> solve(grid(PS_SJ, template)))
    out === nothing && continue
    results[label] = out
    push!(rows, ("Palm Springs–San Jacinto", label, join(size(template), "×"), 12, round(t; digits = 1)))
end
if haskey(results, "900 m")
    out = results["900 m"]
    save_figure("maps_pssj.png", map_figure((
        (hottest(out), "Hottest soil surface (°C)", :thermal, "°C"),
        (coldest(out), "Coldest soil surface (°C)", :thermal, "°C"),
        (deepest_snow(out), "Deepest snow (cm)", Reverse(:ice), "cm"),
        (srtm(PS_SJ, 10), "Elevation (m)", :terrain, "m"),
    ); sites, size = (880, 720)))
end
if haskey(results, "1.8 km") && haskey(results, "900 m")
    save_figure("maps_pssj_resolution.png", map_figure((
        (hottest(results["1.8 km"]), "Hottest soil surface, 1.8 km (°C)", :thermal, "°C"),
        (hottest(results["900 m"]), "Hottest soil surface, 900 m (°C)", :thermal, "°C"),
    ); sites, size = (880, 360)))
end
end

# Mont Aigoual at the native SRTM resolution, the grid of the MicroclimateTalk. The output is cached in
# docs/render/cache/ (not committed) so that the figures can be redrawn without solving again.
run_aigoual = isempty(ARGS) || "aigoual" in ARGS
cache_file = joinpath(mkpath(joinpath(@__DIR__, "cache")), "aigoual.nc")
out = if !run_aigoual
    nothing
elseif isfile(cache_file)
    RasterStack(cache_file)
else
    o, t = timed("Mont Aigoual 66×66", () -> solve(grid(AIGOUAL, SRTM)))
    o === nothing ? nothing : (push!(rows, ("Mont Aigoual", "SRTM, ~90 m", join(size(o.snow_depth)[1:2], "×"), 12, round(t; digits = 1)));
                               write(cache_file, strip_to_canonical(o); force = true); RasterStack(cache_file))
end
if out !== nothing
    dem = srtm(AIGOUAL, 1)
    dem = rebuild(dem; data = Float64.(parent(dem)))
    times = collect(lookup(out, Ti))
    pick(m, h) = findfirst(t -> month(t) == m && hour(t) == h, times)
    # A summer day: morning sun on east-facing slopes, afternoon sun on west-facing ones.
    hours = (8, 11, 14, 17)
    soil = Float64.(out.soil_temperature[depth = 3])        # 5 cm, already °C after strip_to_canonical
    save_figure("maps_aigoual_drape.png", drape_panels(soil, dem, [pick(7, h) for h in hours],
        ["July $(lpad(h, 2, '0')):00" for h in hours]; title = "Soil temperature at 5 cm, Mont Aigoual, July",
        label = "°C", ncols = 2, colorrange = extrema(skipmissing(soil[Ti = [pick(7, h) for h in hours]]))))
    snow = Float64.(out.snow_depth)
    # The months of the melt: the snowfall of the climatology is uniform over this small area, so the pattern
    # comes from the melt, which follows elevation.
    save_figure("maps_aigoual_snow.png", drape_panels(snow, dem, [pick(m, 12) for m in (2, 3, 4, 5)],
        ["February", "March", "April", "May"]; title = "Snow depth at noon, Mont Aigoual", label = "cm",
        colormap = Reverse(:ice), ncols = 2))
end
isempty(rows) || save_table("maps_$(isempty(ARGS) ? "all" : join(ARGS, "_")).csv", ["area", "resolution", "grid", "days", "seconds"], rows)
println("machine: ", MACHINE)
