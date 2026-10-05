# Swapping data sets: one model, weather_source changed, at Madison (2000) and Melbourne (2025).
include(joinpath(@__DIR__, "common.jl"))
const M = MicroclimateMapper

function run_point(source, point, dates; snow = true)
    model = map_model(; micro = micro_model(; snow), weather = source, soil_moisture_source = CPCSoil)
    solve(MicroVectorProblem(; model, points = [point], dates, soil_profile = example_soil_profile(DEPTHS),
        init = start(; snow)))
end

describe(S) = (string(nameof(typeof(M.weather_calendar(S)))),
               M.native_timestep(S) isa M.MinMax ? "daily extremes" : string(M.samples_per_day(M.native_timestep(S)), " a day"))

# Daily maximum and minimum of a layer through its time axis, one value per day solved.
function daily_extremes(series)
    t = collect(lookup(series, Ti)); v = collect(series)
    days = unique(Date.(t))
    (days, [maximum(v[Date.(t) .== d]) for d in days], [minimum(v[Date.(t) .== d]) for d in days])
end

# One untimed run first, so that compilation is not counted in the first timed run.
timed("warm-up (not reported)", () -> run_point(CRUCL2, (-89.40123, 43.07305), Date(2000, 1, 1):Day(1):Date(2000, 12, 31)))

rows = Any[]
for (site, point, dates, sources) in filter(s -> isempty(ARGS) || s[1] in ARGS, (
        ("Madison", (-89.40123, 43.07305), Date(2000, 1, 1):Day(1):Date(2000, 12, 31),
         (CRUCL2, WorldClim{Climate}, TerraClimate{Historical})),
        ("Melbourne", (144.96, -37.81), Date(2025, 1, 1):Day(1):Date(2025, 12, 31),
         (CRUCL2, SILO))))
    fig = Figure(size = (900, 620))
    ax1 = Axis(fig[1, 1]; title = "$site: soil surface, daily maximum", ylabel = "°C")
    ax2 = Axis(fig[2, 1]; title = "$site: air at 1.2 m, daily minimum", ylabel = "°C")
    ax3 = Axis(fig[3, 1]; title = "$site: snow depth, daily maximum", ylabel = "cm")
    colours = Makie.wong_colors()
    for (k, S) in enumerate(sources)
        out, t = timed("$site $S", () -> run_point(S, point, dates))
        cal, step = describe(S)
        out === nothing && (push!(rows, (site, string(S), cal, step, "failed", "")); continue)
        days, tmax, _ = daily_extremes(out.soil_temperature[point = 1, depth = 1])
        _, _, amin = daily_extremes(out.air_temperature[point = 1, height = 2])
        _, smax, _ = daily_extremes(out.snow_depth[point = 1])
        # A monthly data set solves one representative day per month, Microclimate.jl's DEFAULT_DAYS (mid-month).
        x = length(days) <= 12 ? Microclimate.DEFAULT_DAYS[month.(days)] : dayofyear.(days)
        plot! = length(days) <= 12 ? scatterlines! : lines!
        color = colours[k]
        plot!(ax1, x, celsius.(tmax); label = string(S), color)
        plot!(ax2, x, celsius.(amin); color)
        plot!(ax3, x, ustrip.(u"cm", smax); color)
        push!(rows, (site, string(S), cal, step, length(days), round(t; digits = 1)))
    end
    Legend(fig[1:3, 2], ax1; framevisible = false)
    ax3.xlabel = "Day of year"
    site == "Melbourne" && (delete!(ax3); rowsize!(fig.layout, 3, 0); ax2.xlabel = "Day of year")
    save_figure("swap_$(lowercase(site)).png", fig)
end
save_table("swap_data_$(isempty(ARGS) ? "all" : join(ARGS, "_")).csv", ["site", "source", "calendar", "within a day", "days solved", "seconds"], rows)
println("machine: ", MACHINE)
