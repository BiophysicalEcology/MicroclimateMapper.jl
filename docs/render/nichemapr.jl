# NicheMapR's micro_terra at its default site (Madison) for 2000, against the same run here with TerraClimate.
# Reference: test/R/micro_terra_test.R, written to test/data/micro_terra/.
include(joinpath(@__DIR__, "common.jl"))
using DelimitedFiles

refdir = joinpath(@__DIR__, "..", "..", "test", "data", "micro_terra")
read_csv(f) = (d = readdlm(joinpath(refdir, f), ','; header = true); (; data = d[1], header = replace.(vec(d[2]), "\"" => "")))
metout, soil = read_csv("metout_monthly_terra.csv"), read_csv("soil_monthly_terra.csv")
col(t, name) = Float64.(t.data[:, findfirst(==(name), t.header)])

model = map_model(; micro = micro_model(; snow = false, heights = [0.01, 1.2]u"m"), weather = TerraClimate{Historical})
out, t = timed("TerraClimate, Madison, 2000", () -> solve(MicroVectorProblem(; model, points = [(-89.4557, 43.1379)],
    dates = Date(2000, 1, 1):Day(1):Date(2000, 12, 31), soil_profile = example_soil_profile(DEPTHS),
    init = start(; snow = false))))

rmse(a, b) = sqrt(mean((a .- b) .^ 2))
pairs_ = (("air at 1 cm", celsius.(collect(out.air_temperature[point = 1, height = 1])), col(metout, "TALOC")),
          ("air at 1.2 m", celsius.(collect(out.air_temperature[point = 1, height = 2])), col(metout, "TAREF")),
          ("soil surface", celsius.(collect(out.soil_temperature[point = 1, depth = 1])), col(soil, "D0cm")),
          ("soil 10 cm", celsius.(collect(out.soil_temperature[point = 1, depth = 4])), col(soil, "D10cm")),
          ("soil 50 cm", celsius.(collect(out.soil_temperature[point = 1, depth = 8])), col(soil, "D50cm")))
save_table("nichemapr_terra.csv", ["variable", "RMSE (K)", "mean difference (K)"],
    [(name, round(rmse(here, ref); digits = 2), round(mean(here .- ref); digits = 2)) for (name, here, ref) in pairs_])

fig = Figure(size = (900, 520))
for (k, (name, here, ref)) in enumerate(pairs_[[1, 3, 4, 5]])
    row, c = fldmod1(k, 2)
    ax = Axis(fig[row, c]; title = name, ylabel = c == 1 ? "°C" : "", xlabel = row == 2 ? "Hour of the 12 representative days" : "")
    lines!(ax, ref; color = :grey60, linewidth = 3, label = "NicheMapR micro_terra")
    lines!(ax, here; color = :darkorange, label = "MicroclimateMapper.jl")
end
Legend(fig[3, 1:2], content(fig[1, 1]); orientation = :horizontal, framevisible = false)
save_figure("nichemapr_terra.png", fig)
println("solve ", round(t; digits = 1), " s; machine: ", MACHINE)
