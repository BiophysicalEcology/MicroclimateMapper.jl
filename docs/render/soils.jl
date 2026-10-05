# Soils: a profile from SLGA texture at Melbourne, three pedotransfer functions, and soil moisture solved from SILO rainfall.
include(joinpath(@__DIR__, "common.jl"))

const YEAR = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)
const MELBOURNE = (144.96, -37.81)
d = ustrip.(u"cm", DEPTHS)

ptfs = (CosbyUnivariate(), CosbyMultivariate(), Campbell1985())
profiles = []
rows = Any[]
for ptf in ptfs
    p, _ = timed("SLGA $(nameof(typeof(ptf)))", () -> build_soil_profile(SLGA, MELBOURNE; depths = DEPTHS, pedotransfer_model = ptf))
    p === nothing && continue
    push!(profiles, (string(nameof(typeof(ptf))), p))
    push!(rows, (nameof(typeof(ptf)), round(p.campbell_b[1]; digits = 2), round(ustrip(u"J/kg", p.air_entry_potential[1]); digits = 2),
                 round(ustrip(u"kg*s/m^3", p.saturated_conductivity[1]); sigdigits = 3),
                 round(ustrip(u"m^3/m^3", p.field_capacity[1]); digits = 3), round(ustrip(u"m^3/m^3", p.wilting_point[1]); digits = 3)))
end
save_table("soils_surface.csv", ["pedotransfer", "Campbell b", "air-entry potential (J/kg)", "saturated conductivity (kg s/m³)",
    "field capacity (m³/m³)", "wilting point (m³/m³)"], rows)

fig = Figure(size = (900, 380))
for (k, (label, f)) in enumerate((("Campbell b", p -> p.campbell_b), ("Air-entry potential (J/kg)", p -> ustrip.(u"J/kg", p.air_entry_potential)),
                                  ("Saturated conductivity (kg s m⁻³)", p -> ustrip.(u"kg*s/m^3", p.saturated_conductivity))))
    ax = Axis(fig[1, k]; xlabel = label, ylabel = k == 1 ? "Depth (cm)" : "", yreversed = true)
    for (name, p) in profiles
        scatterlines!(ax, f(p), d; label = name)
    end
    k == 1 && axislegend(ax; position = :rb, labelsize = 9)
end
save_figure("soils_profiles.png", fig)

# Soil moisture solved from SILO's daily rainfall at Melbourne for 2025, with Microclimate.jl's example loam and the
# SLGA profile under two pedotransfer functions. With prescribed moisture the hydraulic properties do not enter;
# solved, they set how fast rain soaks in and drains.
# Not yet on the page (October 2026): with the SLGA profiles, Microclimate.jl's dynamic soil moisture dries the top
# layers to 0% and does not rewet them, which is physically impossible; the example loam behaves. Report, then use.
const YEAR_2025 = Date(2025, 1, 1):Day(1):Date(2025, 12, 31)
dynamic = map_model(; weather = SILO,
    micro = micro_model(; snow = false, config = MicroConfig(; soil_moisture_strategy = DynamicSoilMoisture())))
soils = [("example loam", example_soil_profile(DEPTHS))]
for (name, p) in profiles
    name in ("CosbyMultivariate", "Campbell1985") && push!(soils, ("SLGA, $name", p.soil_profile))
end
runs = []
for (label, sp) in soils
    out, t = timed("SILO 2025, dynamic soil moisture, $label", () -> solve(MicroVectorProblem(; model = dynamic,
        points = [MELBOURNE], dates = YEAR_2025, soil_profile = sp, init = start(; snow = false))))
    out === nothing && continue
    push!(runs, (label, out))
end
if !isempty(runs)
    fig = Figure(size = (900, 700))
    daily(series, f) = (t = collect(lookup(series, Ti)); v = collect(series); days = unique(Date.(t));
                        (dayofyear.(days), [f(v[Date.(t) .== d]) for d in days]))
    panels = (("Soil moisture at 5 cm, daily mean", "% by volume", o -> o.soil_moisture[point = 1, depth = 3], v -> 100mean(v)),
              ("Soil moisture at 50 cm, daily mean", "% by volume", o -> o.soil_moisture[point = 1, depth = 8], v -> 100mean(v)),
              ("Soil surface, daily maximum", "°C", o -> o.soil_temperature[point = 1, depth = 1], v -> ustrip(u"°C", maximum(v))))
    axes = [Axis(fig[k, 1]; title, ylabel, xlabel = k == 3 ? "Day of year, 2025" : "") for (k, (title, ylabel, _, _)) in enumerate(panels)]
    for (j, (label, out)) in enumerate(runs), (ax, (_, _, layer, f)) in zip(axes, panels)
        x, y = daily(layer(out), f)
        lines!(ax, x, y; label, color = Makie.wong_colors()[j])
    end
    linkxaxes!(axes...)
    Legend(fig[1:3, 2], axes[1]; framevisible = false)
    save_figure("soils_moisture.png", fig)
    summary = Any[]
    for (label, out) in runs
        m5 = 100 .* collect(out.soil_moisture[point = 1, depth = 3]); m50 = 100 .* collect(out.soil_moisture[point = 1, depth = 8])
        ts = celsius.(collect(out.soil_temperature[point = 1, depth = 1]))
        push!(summary, (label, round(mean(m5); digits = 1), round(minimum(m5); digits = 1), round(maximum(m5); digits = 1),
                        round(mean(m50); digits = 1), round(maximum(ts); digits = 1)))
    end
    save_table("soils_moisture.csv", ["soil", "mean at 5 cm (%)", "min at 5 cm (%)", "max at 5 cm (%)", "mean at 50 cm (%)",
        "hottest surface (°C)"], summary)
    foreach(println, summary)
end
println("machine: ", MACHINE)
