# Canopy on Saba: canopy heights from the LiDAR surface model less the terrain model, as in SolarRadiation.jl's
# Vegetation section; one run per canopy-height class, with the multilayer canopy of Microclimate.jl.
include(joinpath(@__DIR__, "common.jl"))
using Downloads
using Microclimate: DEFAULT_DEPTHS

const YEAR = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)
data_dir = mkpath(joinpath(ENV["RASTERDATASOURCES_PATH"], "saba"))
for f in ("saba.tif", "saba_dsm.tif")
    isfile(joinpath(data_dir, f)) ||
        Downloads.download("https://github.com/Deltares/Geomorphometry.jl/releases/download/v0.6.0/$f", joinpath(data_dir, f))
end

# Canopy height of every land cell, clamped to 0–30 m as in SolarRadiation.jl's tutorial.
dtm = Raster(joinpath(data_dir, "saba.tif"))
dsm = Raster(joinpath(data_dir, "saba_dsm.tif"))
canopy = [ismissing(g) || ismissing(s) ? missing : clamp(Float64(s - g), 0.0, 30.0) for (g, s) in zip(parent(dtm), parent(dsm))]
heights_on_land = collect(skipmissing(canopy))
quantiles = quantile(heights_on_land, [0.25, 0.5, 0.75, 0.95])
println("canopy height quantiles (m): ", round.(quantiles; digits = 1))

# Height grid and plant-area profile, after MicroclimateTests.jl/demos/monthly_canopy.jl.
function canopy_heights(canopy_height, reference_height)
    near = [0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.0, 1.5, 2.0]
    h, r = ustrip(u"m", canopy_height), ustrip(u"m", reference_height)
    inside = sort(unique(vcat(filter(<=(h), near), range(0.0, h; length = 11)[2:end])))
    above = collect(range(h, r; length = 7)[2:end])
    sort(unique(vcat(inside, above))) .* u"m"
end
function pai_profile(heights, canopy_height, total)
    n = count(<=(canopy_height), heights)
    (; layer_heights) = Microclimate.canopy_layer_heights(heights, canopy_height, n)
    raw = fill(1.0 / u"m", length(layer_heights))
    raw .* (total / sum(raw))   # as plant_area_index_profile in MicroclimateTests.jl/demos/monthly_canopy.jl
end

function canopy_model(canopy_height, plant_area_index_total)
    reference_height = max(2.0u"m", canopy_height + 2.0u"m")
    heights = canopy_heights(canopy_height, reference_height)
    plant_area_index = pai_profile(heights, canopy_height, plant_area_index_total)
    model = MultilayerCanopy(; canopy_height, plant_area_index, shortwave_model = TwoStreamRadiation(),
        convergence_model = PicardCanopyConvergence(;
            convergence = IterationToleranceConvergence(; tolerance = 0.1u"K", max_iterations_per_day = 20), relaxation = 0.7))
    micro = MicroModel(; hours = 0:1:23, depths = DEPTHS, heights,
        soil_properties_model = example_soil_properties_model(), soil_hydraulic_model = example_soil_hydraulic_model(),
        snow_model = NoSnow(), canopy_model = model,
        config = MicroConfig(; convergence = IterationToleranceConvergence(; tolerance = 0.1u"K", max_iterations_per_day = 20),
            soil_moisture_strategy = PrescribedSoilMoisture(),
            canopy_soil_convergence = IterationToleranceConvergence(; tolerance = 0.1u"K", max_iterations_per_day = 20)))
    map_model(; micro, weather = TerraClimate{Historical}, output_layers = (
        LayerSpec(:soil_temperature, :soil), LayerSpec(:air_temperature, :profile),
        LayerSpec(:leaf_temperature, :canopy), LayerSpec(:canopy_air_temperature, :canopy, :air_temperature)))
end
open_model = map_model(; micro = micro_model(; snow = false), weather = TerraClimate{Historical})

point = [(-63.233, 17.635)]   # the upper slopes of Mount Scenery
classes = [("open ground", 0.0), ("forest", round(quantiles[4]; digits = 0))]
results = []
rows = Any[]
for (name, h) in classes
    model = h == 0 ? open_model : canopy_model(h * u"m", h < 5 ? 1.5 : 4.0)
    out, t = timed("$name, $(h) m", () -> solve(MicroVectorProblem(; model, points = point, dates = YEAR,
        soil_profile = example_soil_profile(DEPTHS), init = start(; snow = false))))
    out === nothing && continue
    push!(results, (name, h, out))
    push!(rows, (name, h, round(t; digits = 1)))
end

fig = Figure(size = (900, 360))
ax1 = Axis(fig[1, 1]; title = "Soil surface, June", xlabel = "Hour", ylabel = "°C")
ax2 = Axis(fig[1, 2]; title = "Soil at 10 cm, June", xlabel = "Hour")
for (name, h, out) in results
    june = findall(t -> month(t) == 6, collect(lookup(out, Ti)))
    lines!(ax1, 0:23, celsius.(collect(out.soil_temperature[point = 1, depth = 1])[june]); label = h == 0 ? name : "$name, $(Int(h)) m")
    lines!(ax2, 0:23, celsius.(collect(out.soil_temperature[point = 1, depth = 4])[june]))
end
Legend(fig[1, 3], ax1; framevisible = false)
save_figure("canopy_saba.png", fig)

# Leaves of the forest: leaf temperature in the top, middle and bottom canopy layers (layer 1 is the top), each
# against the air at the same height, on the representative day of June.
forest = findfirst(r -> r[2] > 0, results)
if forest !== nothing
    _, h, out = results[forest]
    n = size(out.leaf_temperature, Dim{:canopy_layer})
    layer_heights = ustrip.(u"m", Microclimate.canopy_layer_heights(
        canopy_heights(h * u"m", max(2.0u"m", h * u"m" + 2.0u"m")), h * u"m", n).layer_heights)
    june = findall(t -> month(t) == 6, collect(lookup(out, Ti)))
    leaves = Figure(size = (900, 360))
    ax = Axis(leaves[1, 1]; title = "Forest, $(Int(h)) m, June: leaves (solid) and air (dashed)", xlabel = "Hour", ylabel = "°C")
    for (k, (layer, name)) in enumerate(((1, "top"), (cld(n, 2), "middle"), (n, "bottom")))
        color = Makie.wong_colors()[k]
        label = "$name layer, $(round(layer_heights[layer]; sigdigits = 2)) m"
        lines!(ax, 0:23, celsius.(collect(out.leaf_temperature[point = 1, canopy_layer = layer])[june]); color, label)
        lines!(ax, 0:23, celsius.(collect(out.canopy_air_temperature[point = 1, canopy_layer = layer])[june]); color, linestyle = :dash)
    end
    ax2 = Axis(leaves[1, 2]; title = "Leaf less air", xlabel = "Hour", ylabel = "K")
    for (k, layer) in enumerate((1, cld(n, 2), n))
        Δ = ustrip.(u"K", collect(out.leaf_temperature[point = 1, canopy_layer = layer])[june]) .-
            ustrip.(u"K", collect(out.canopy_air_temperature[point = 1, canopy_layer = layer])[june])
        lines!(ax2, 0:23, Δ; color = Makie.wong_colors()[k])
    end
    hlines!(ax2, 0.0; color = :grey60)
    Legend(leaves[1, 3], ax; framevisible = false)
    save_figure("canopy_saba_leaves.png", leaves)
    top = ustrip.(u"K", collect(out.leaf_temperature[point = 1, canopy_layer = 1])[june]) .-
          ustrip.(u"K", collect(out.canopy_air_temperature[point = 1, canopy_layer = 1])[june])
    bottom = ustrip.(u"K", collect(out.leaf_temperature[point = 1, canopy_layer = n])[june]) .-
             ustrip.(u"K", collect(out.canopy_air_temperature[point = 1, canopy_layer = n])[june])
    println("leaf less air, top layer: max ", round(maximum(top); digits = 2), " K at ", argmax(top) - 1, ":00, min ",
        round(minimum(top); digits = 2), " K; bottom layer: max ", round(maximum(bottom); digits = 2), ", min ", round(minimum(bottom); digits = 2))
end

hist = Figure(size = (520, 320))
ax = Axis(hist[1, 1]; xlabel = "Canopy height (m)", ylabel = "Land cells", title = "Saba: surface model less terrain model")
CairoMakie.hist!(ax, heights_on_land; bins = 0:1:30)
save_figure("canopy_saba_heights.png", hist)
save_table("canopy.csv", ["class", "canopy height (m)", "seconds"], rows)
println("machine: ", MACHINE)
