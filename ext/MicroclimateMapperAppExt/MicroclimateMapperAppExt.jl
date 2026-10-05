module MicroclimateMapperAppExt

# A local interactive app for point microclimates: pick a place on a map, swap data sets and settings, compare data
# sets, and copy a script that reproduces the run. `using MicroclimateMapper, Bonito, WGLMakie; app()`.

using MicroclimateMapper
using MicroclimateMapper: geocode
using Microclimate
using Microclimate: example_soil_profile, example_soil_properties_model, example_soil_hydraulic_model
using RasterDataSources
import Rasters
using Rasters: Raster, Ti, X, Y, crop, lookup, rebuild
using Bonito
using WGLMakie
using Dates
using Unitful
using SolarRadiation: SolarRadiation, solar_geometry
import MicroclimateMapper: app

const LEAFLET_JS = Asset("https://unpkg.com/leaflet@1.9.4/dist/leaflet.js")
const LEAFLET_CSS = Asset("https://unpkg.com/leaflet@1.9.4/dist/leaflet.css")

# ── Data sets offered ─────────────────────────────────────────────────────────────────────────────────────────────

struct Choice
    label::String
    source::Any
    code::String                                 # how the source is written in a script
    region::NTuple{4,Float64}                    # lon min, lon max, lat min, lat max
    years::Union{Nothing,UnitRange{Int}}         # nothing: a climatology
    note::String
end

const GLOBAL = (-180.0, 180.0, -90.0, 90.0)
const USA = (-125.0, -66.5, 24.5, 49.5)
const AUSTRALIA = (112.0, 154.0, -44.0, -10.0)
const CHOICES = [
    Choice("CRU CL 2.0 (1961–1990 climatology)", CRUCL2, "CRUCL2", GLOBAL, nothing,
        "monthly normals, 18 km; one small download; land only"),
    Choice("WorldClim (1970–2000 climatology)", WorldClim{Climate}, "WorldClim{Climate}", GLOBAL, nothing,
        "monthly normals, 18 km; land only; not corrected for elevation"),
    Choice("TerraClimate (a year, monthly)", TerraClimate{Historical}, "TerraClimate{Historical}", GLOBAL, 1958:2024,
        "monthly, 4 km; large yearly downloads; not corrected for elevation"),
    Choice("gridMET (a year, daily; USA)", GRIDMET, "GRIDMET", USA, 1979:2025,
        "daily, 4 km, contiguous USA; large yearly downloads; not corrected for elevation"),
    Choice("SILO (a year, daily; Australia)", SILO, "SILO", AUSTRALIA, 1889:2025,
        "daily, 5 km, Australia; needs ENV[\"SILO_EMAIL\"]; not corrected for elevation"),
    Choice("NCEP reanalysis (a year, 6-hourly)", NCEP{SurfaceFlux, 1}, "NCEP{SurfaceFlux, 1}", GLOBAL, 1948:2025,
        "6-hourly, ~1.9°; not corrected for elevation"),
    Choice("ERA5 reanalysis (a year, hourly)", ERA5, "ERA5", GLOBAL, 1940:2025,
        "hourly, 0.25°; needs `using ZarrDatasets`; slow for a year without a CDS key"),
    Choice("ERA5-Land (a year, hourly)", ERA5ECMWFLand, "ERA5ECMWFLand", GLOBAL, 1950:2025,
        "hourly, 0.1° (~9 km), land only, cloud from ERA5; read from ECMWF's store: needs a CDS key and `using ZarrDatasets`"),
]

# ERA5 by its fast route: ECMWF's geo-chunked store, read point by point through PointDataSources.jl. It needs a
# Copernicus (CDS) API key and `using ZarrDatasets`; without a key, ERA5 is read from Google's ARCO store, slowly.
# ERA5-Land is read from ECMWF's gridded stores by RasterDataSources.jl, which takes the same key under another name.
function use_cds_era5!(key::AbstractString)
    isempty(strip(key)) && return false
    ENV["CDS_API_KEY"] = strip(key)    # PointDataSources.jl
    ENV["CDSAPI_KEY"] = strip(key)     # RasterDataSources.jl
    @eval MicroclimateMapper.loader(::Type{<:ERA5}) = MicroclimateMapper.PointQuery()
    return true
end

covers(c::Choice, lon, lat) = c.region[1] <= lon <= c.region[2] && c.region[3] <= lat <= c.region[4]

# ── Settings and the run ──────────────────────────────────────────────────────────────────────────────────────────

# Every run has these depths and heights; the app chooses which to plot. The top height is the reference height.
const DEPTHS = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
const HEIGHTS = [0.01, 0.05, 0.10, 0.30, 0.50, 1.0, 2.0]u"m"
length_label(x) = x < 1u"m" ? string(round(Int, ustrip(u"cm", x)), " cm") : string(round(Int, ustrip(u"m", x)), " m")
const DEPTH_LABELS = [x == 0u"cm" ? "surface" : length_label(x) for x in DEPTHS]
const HEIGHT_LABELS = length_label.(HEIGHTS)

Base.@kwdef struct Settings
    lon::Float64 = -116.545
    lat::Float64 = 33.830
    year::Int = 2010
    shade::Float64 = 0.0           # fraction
    snow::Bool = false
    dynamic_moisture::Bool = false # soil moisture solved from rainfall; otherwise prescribed
    terrain::Bool = false          # slope, aspect and horizon from SRTM; otherwise the three below, at CRUCL2's elevation
    slope::Float64 = 0.0           # degrees
    aspect::Float64 = 180.0        # degrees clockwise from north
    horizon::Float64 = 0.0         # degrees, the same in every direction
    albedo::Float64 = 0.15
    soil::Symbol = :example        # :example, :SoilGrids, :SLGA
    data_path::String = get(ENV, "RASTERDATASOURCES_PATH", "")
end

dates_for(c::Choice, s::Settings) = c.years === nothing ? (Date(2000, 1, 1):Day(1):Date(2000, 12, 31)) :
    (Date(s.year, 1, 1):Day(1):Date(s.year, 12, 31))

function soil_profile(s::Settings)
    s.soil == :example && return example_soil_profile(DEPTHS)
    source = s.soil == :SLGA ? SLGA : SoilGrids
    build_soil_profile(source, (s.lon, s.lat); depths = DEPTHS).soil_profile
end

# Two caches make repeated runs fast. Terrain (elevation, slope, aspect, horizons) depends only on the DEM and the
# point, so it is reused whatever else changes; it is most of `init` after the first run. A result is reused when the
# same data set is run again with the same settings, as when another data set is added to a comparison.
const TERRAIN = Dict{Tuple{Float64,Float64,Bool},Any}()
const RESULTS = Dict{Tuple{String,Settings},Any}()
const RUN_LOCK = ReentrantLock()     # one Run at a time
const CACHE_LOCK = ReentrantLock()   # the caches are shared by the data sets solved in parallel

function solve_point(c::Choice, s::Settings)
    hit = lock(() -> get(RESULTS, (c.code, s), nothing), CACHE_LOCK)
    hit === nothing || return hit
    micro_model = MicroModel(; depths = DEPTHS, heights = HEIGHTS,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model = example_soil_hydraulic_model(),
        snow_model = s.snow ? SnowModel() : NoSnow(),
        config = MicroConfig(; soil_moisture_strategy = s.dynamic_moisture ? DynamicSoilMoisture() : PrescribedSoilMoisture()))
    model = MicroMapModel(; micro_model, dem_source = dem_source(s), weather_source = c.source,
        soil_moisture_source = c.source === TerraClimate{Historical} ? nothing : CPCSoil,
        surface_albedo_source = s.albedo, roughness_height_source = 0.004u"m",
        compute_terrain = s.terrain,
        output_layers = OUTPUT_LAYERS)
    start = s.snow ? (; soil_moisture = fill(0.2, length(DEPTHS)), snow_depth = 0.0u"cm") :
                     (; soil_moisture = fill(0.2, length(DEPTHS)))
    problem(data) = MicroVectorProblem(; model, points = [(s.lon, s.lat)], dates = dates_for(c, s),
        soil_profile = soil_profile(s), init = start, data)
    key = (s.lon, s.lat, s.terrain)
    known = lock(() -> get(TERRAIN, key, nothing), CACHE_LOCK)
    cache = nothing
    if known === nothing   # first run here: load the terrain once (CRUCL2's elevation and flat, or all of it from SRTM)
        first_cache = init(problem((; shade = s.shade)))
        known = terrain(first_cache)
        lock(() -> (TERRAIN[key] = deepcopy(known)), CACHE_LOCK)
        # With terrain from SRTM this cache is already the run; only the sliders' terrain needs a second init.
        s.terrain && (cache = first_cache)
    end
    if cache === nothing
        # Slope, aspect and horizon from the sliders must be in the terrain when `init` runs: init prepares the point.
        t = s.terrain ? deepcopy(known) : set_terrain!(deepcopy(known), s)
        cache = init(problem((; shade = s.shade, terrain = t)))
    end
    out = solve!(cache)
    lock(() -> (RESULTS[(c.code, s)] = out), CACHE_LOCK)
    return out
end

# Off, the elevation comes from CRUCL2's 10-minute grid, which loads in a moment; on, everything comes from SRTM.
dem_source(s::Settings) = s.terrain ? SRTM : CRUCL2

const OUTPUT_LAYERS = (MicroclimateMapper._DEFAULT_OUTPUT_LAYERS...,
    LayerSpec(:soil_water_potential, :soil), LayerSpec(:soil_humidity, :soil), LayerSpec(:diffuse_fraction, :scalar))

# Slope, aspect and horizon set by hand in a copy of the point's terrain, which is then passed to `init` (elevation
# stays that of the DEM at the point, for the elevation correction of the weather).
function set_terrain!(t, s::Settings)
    t.slope .= s.slope * u"°"
    t.aspect .= s.aspect * u"°"
    n = length(first(t.horizon_angles))
    t.horizon_angles .= Ref(MicroclimateMapper.SVector{n}(fill(s.horizon * u"°", n)))
    return t
end

# ── The script that reproduces a run ──────────────────────────────────────────────────────────────────────────────

# The lines of the script that set the terrain, indented to sit inside `run_with`.
function terrain_lines(s::Settings)
    s.terrain && return "# slope, aspect and horizon from the SRTM DEM (about 90 m)"
    lines = [
        "# Flat terrain at the elevation of CRUCL2's 10-minute cell; slope, aspect and horizon set by hand and passed back in:",
        "t = terrain(init(problem))",
        "t.slope .= $(s.slope)u\"°\"",
        "t.aspect .= $(s.aspect)u\"°\"   # clockwise from north",
        "t.horizon_angles .= Ref(MicroclimateMapper.SVector{32}(fill($(s.horizon)u\"°\", 32)))   # in every direction",
        "problem = MicroVectorProblem(; problem.model, problem.points, problem.dates, problem.soil_profile,",
        "    problem.init, data = (; problem.data..., terrain = t))",
    ]
    return join(lines, "
    ")
end

const SILO_LINES = "# SILO needs ENV[\"SILO_EMAIL\"] set to your email address\n"
const ERA5_LAND_LINES = """using ZarrDatasets   # ERA5-Land, read from ECMWF's store
# ERA5-Land needs your Copernicus (CDS) API key in ENV["CDSAPI_KEY"] or a ~/.cdsapirc file
"""
const ERA5_LINES = """using ZarrDatasets   # ERA5
# With a Copernicus key in ENV["CDS_API_KEY"], read ERA5 point by point from ECMWF's store (fast):
MicroclimateMapper.loader(::Type{<:ERA5}) = MicroclimateMapper.PointQuery()
"""

function script(choices::Vector{Choice}, s::Settings)
    soil = s.soil == :example ? "example_soil_profile(depths)" :
        "build_soil_profile($(s.soil), point; depths).soil_profile"
    moisture(x) = x.source === TerraClimate{Historical} ? "nothing" : "CPCSoil"
    dates(x) = x.years === nothing ? "Date(2000, 1, 1):Day(1):Date(2000, 12, 31)" :
        "Date($(s.year), 1, 1):Day(1):Date($(s.year), 12, 31)"
    """
    # Made by MicroclimateMapper.app() on $(Dates.format(now(), "yyyy-mm-dd HH:MM")).
    using MicroclimateMapper, Microclimate, RasterDataSources, Rasters, Dates, Unitful
    using Microclimate: example_soil_profile, example_soil_properties_model, example_soil_hydraulic_model
    $(any(x -> x.source === ERA5, choices) ? ERA5_LINES : "")$(any(x -> x.source === ERA5ECMWFLand, choices) ? ERA5_LAND_LINES : "")$(any(x -> x.source === SILO, choices) ? SILO_LINES : "")
    ENV["RASTERDATASOURCES_PATH"] = $(repr(s.data_path))   # where data sets are downloaded and cached

    point = ($(round(s.lon; digits = 4)), $(round(s.lat; digits = 4)))   # longitude, latitude
    depths = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
    heights = [0.01, 0.05, 0.10, 0.30, 0.50, 1.0, 2.0]u"m"   # the top one is the reference height
    micro_model = MicroModel(; depths, heights,
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model = example_soil_hydraulic_model(),
        snow_model = $(s.snow ? "SnowModel()" : "NoSnow()"),
        config = MicroConfig(; soil_moisture_strategy = $(s.dynamic_moisture ? "DynamicSoilMoisture()" : "PrescribedSoilMoisture()")))

    # A climatology (CRUCL2, WorldClim) solves 12 representative days whatever year its dates are in.
    function run_with(source, dates; soil_moisture_source = CPCSoil)
        problem = MicroVectorProblem(;
            model = MicroMapModel(; micro_model, dem_source = $(s.terrain ? "SRTM" : "CRUCL2"), weather_source = source, soil_moisture_source,
                surface_albedo_source = $(s.albedo), roughness_height_source = 0.004u"m", compute_terrain = $(s.terrain)),
            points = [point], dates, soil_profile = $soil,
            init = $(s.snow ? "(; soil_moisture = fill(0.2, length(depths)), snow_depth = 0.0u\"cm\")" : "(; soil_moisture = fill(0.2, length(depths)))"),
            data = (; shade = $(s.shade)))
        $(terrain_lines(s))
        cache = init(problem)
        return solve!(cache)
    end

    outputs = [
    $(join(("    run_with($(x.code), $(dates(x)); soil_moisture_source = $(moisture(x)))," for x in choices), "\n"))
    ]

    # Each output is a RasterStack, selected by name: soil temperature at 10 cm, air temperature at 30 cm
    outputs[1].soil_temperature[point = 1, depth = Near(0.1)]
    outputs[1].air_temperature[point = 1, height = Near(0.3)]
    """
end

# ── Figures ───────────────────────────────────────────────────────────────────────────────────────────────────────

const MONTH_NAMES = ["January", "February", "March", "April", "May", "June", "July", "August", "September",
    "October", "November", "December"]
const MONTH_ABBREVIATIONS = first.(MONTH_NAMES, 3)

celsius(x) = ustrip(u"°C", x)
percent(x) = (v = Float64(ustrip(x)); v <= 1.0 ? 100v : v)   # relative humidity may come as a fraction
suction(x) = max(-ustrip(u"J/kg", x), 1e-3)                   # −ψ, positive, for a log scale
kg_m2(x) = ustrip(u"kg/m^2", x)

is_climatology(o) = length(unique(Date.(lookup(o, Ti)))) <= 12

# What part of the run to plot: the whole year, a month or a day.
struct View
    period::Symbol   # :year, :month or :day
    month::Int
    day::Int
end

# The hours to show from one output and their x positions, with the x-axis label. A climatology's 12 representative
# days are joined end to end for the year; its month and its day are both that month's representative day. When the
# first data set is a climatology, the others are shown on the same days.
function hours_to_show(o, v::View, reference_dates)
    t = collect(lookup(o, Ti))
    if reference_dates !== nothing
        if v.period == :year
            i = findall(x -> Date(x) in reference_dates, t)
            return i, Float64[(findfirst(==(Date(t[k])), reference_dates) - 1) * 24 + hour(t[k]) for k in i],
                "The 12 representative days, hour by hour"
        end
        i = findall(x -> Date(x) == reference_dates[v.month], t)
        return i, Float64.(hour.(t[i])), "Hour of the representative day of $(MONTH_NAMES[v.month])"
    end
    v.period == :year && return collect(eachindex(t)), dayofyear.(t) .+ hour.(t) ./ 24, "Day of year"
    if v.period == :month
        i = findall(x -> month(x) == v.month, t)
        return i, day.(t[i]) .+ hour.(t[i]) ./ 24, "Day of $(MONTH_NAMES[v.month])"
    end
    i = findall(x -> month(x) == v.month && day(x) == v.day, t)
    return i, Float64.(hour.(t[i])), "Hour of $(v.day) $(MONTH_NAMES[v.month])"
end

# The day of year a day view shows, or nothing for a year or month of daily data.
function day_shown(v::View, climatology::Bool, year)
    climatology && v.period != :year && return Microclimate.DEFAULT_DAYS[v.month]   # the representative day solved
    v.period == :day && return dayofyear(Date(year, v.month, min(v.day, daysinmonth(Date(year, v.month)))))
    return nothing
end

# Sunrise and sunset in solar time, as the runs use, from SolarRadiation.jl: half the day length is the sunrise hour
# angle (McCullough and Porter 1971, eq. 7), so the sun rises at 12 − H and sets at 12 + H.
declination(lat, doy) = solar_geometry(SolarRadiation.McCulloughPorterSolarGeometry(), lat * u"°";
    day_of_year = doy, hour_angle = 0.0u"rad").solar_declination
function sun_hours(lat, doy)
    H = SolarRadiation.sunrise_hour_angle(declination(lat, doy), lat * u"°").H₋
    return 12 - H, 12 + H
end

# The sun's path through a day: azimuth (clockwise from north) and elevation above the horizon, in degrees.
function sun_path(lat, doy)
    δ = declination(lat, doy)
    path = Tuple{Float64,Float64}[]
    for t in 0:0.1:24
        h = (π / 12) * (t - 12) * u"rad"
        z = solar_geometry(SolarRadiation.McCulloughPorterSolarGeometry(), lat * u"°"; day_of_year = doy, hour_angle = h).zenith_angle
        elevation = 90 - rad2deg(ustrip(u"rad", z))
        # solar_azimuth_angle returns radians as a plain number, or degrees at solar noon: convert both to degrees.
        azimuth = rad2deg(Float64(uconvert(NoUnits, SolarRadiation.solar_azimuth_angle(h, lat * u"°", δ))))
        elevation > 0 && push!(path, (mod(azimuth, 360), elevation))
    end
    return path
end

# Plot groups. Each draws one or more panels; soil water and snow are offered only when their model option is on.
const GROUPS = [
    (:soil_temperature, "Soil temperature"),
    (:above_ground, "Above ground: air temperature, wind speed, humidity"),
    (:radiation, "Radiation: solar, longwave down and up"),
    (:soil_water, "Soil water: moisture, humidity, water potential"),
    (:snow, "Snow"),
    (:dew_frost, "Dew and frost at the ground surface"),
]
const DEFAULT_GROUPS = (:soil_temperature, :above_ground)
const STEFAN_BOLTZMANN = 5.670374419e-8   # W m⁻² K⁻⁴
const SURFACE_EMISSIVITY = 0.98           # as Microclimate.jl's default

available(group, s::Settings) = group == :soil_water ? s.dynamic_moisture : group == :snow ? s.snow : true

# Solar radiation on the slope, as the soil surface receives it. `global_radiation` is on a horizontal surface (the
# horizon has already removed the direct beam when the sun is behind it). The soil energy balance of Microclimate.jl
# scales the direct beam by cos(zenith on the slope) / cos(zenith) and the diffuse part by the sky view factor; the
# same is done here hour by hour with the sun's position from SolarRadiation.jl. `nothing` on flat ground.
function slope_solar(out, s::Settings)
    t = point_terrain(s)
    t === nothing && return nothing
    slope = ustrip(u"°", first(t.slope))
    slope > 0 || return nothing
    aspect = ustrip(u"°", first(t.aspect))
    sky_view = MicroclimateMapper._sky_view_from_horizon(first(t.horizon_angles))
    global_radiation = ustrip.(u"W/m^2", collect(out.global_radiation[point = 1]))
    diffuse = collect(out.diffuse_fraction[point = 1])
    times = collect(lookup(out, Ti))
    return map(eachindex(global_radiation)) do i
        h = (π / 12) * (hour(times[i]) + minute(times[i]) / 60 - 12) * u"rad"
        geometry = solar_geometry(SolarRadiation.McCulloughPorterSolarGeometry(), s.lat * u"°";
            day_of_year = dayofyear(times[i]), hour_angle = h)
        zenith = rad2deg(ustrip(u"rad", geometry.zenith_angle))
        zenith >= 90 && return global_radiation[i]
        azimuth = rad2deg(Float64(uconvert(NoUnits, SolarRadiation.solar_azimuth_angle(h, s.lat * u"°", geometry.solar_declination))))
        cos_slope_zenith = cosd(zenith) * cosd(slope) + sind(zenith) * sind(slope) * cosd(azimuth - aspect)
        direct = global_radiation[i] * (1 - diffuse[i])
        max(0.0, direct / cosd(zenith) * cos_slope_zenith) + global_radiation[i] * diffuse[i] * sky_view
    end
end

function figure(results, v::View, depths_shown, heights_shown, groups, s::Settings)
    c, out = first(results)
    climatology = is_climatology(out)
    reference_dates = climatology ? unique(Date.(lookup(out, Ti))) : nothing
    idx, x, xlabel = hours_to_show(out, v, reference_dates)
    isempty(idx) && return DOM.p("Nothing to show for that day.")
    compare = length(results) > 1
    h1 = isempty(heights_shown) ? 1 : first(heights_shown)
    groups = [g for g in groups if available(g, s)]
    # Each panel: title, y label, a function drawing into an axis, the y scale. Legends go beside the panels.
    panels = Tuple{String,String,Function,Any}[]
    soil_panel!(title, ylabel, layer, f; yscale = identity) = push!(panels, ("$title, $(c.label)", ylabel, ax -> begin
        for d in depths_shown
            lines!(ax, x, f.(collect(getproperty(out, layer)[point = 1, depth = d]))[idx]; label = DEPTH_LABELS[d])
        end
    end, yscale))
    air_panel!(title, ylabel, layer, f) = push!(panels, (compare ? "$title at $(HEIGHT_LABELS[h1])" : "$title, $(c.label)", ylabel, ax -> begin
        if compare
            for (k, (ck, o)) in enumerate(results)
                i, xk, _ = hours_to_show(o, v, reference_dates)
                isempty(i) || lines!(ax, xk, f.(collect(getproperty(o, layer)[point = 1, height = h1]))[i];
                    color = Makie.wong_colors()[mod1(k, 7)], label = ck.label)
            end
        else
            for h in heights_shown
                lines!(ax, x, f.(collect(getproperty(out, layer)[point = 1, height = h]))[idx]; label = HEIGHT_LABELS[h])
            end
        end
    end, identity))
    function series!(ax, o, layer, f; kw...)
        i, xk, _ = hours_to_show(o, v, reference_dates)
        isempty(i) || lines!(ax, xk, f.(collect(getproperty(o, layer)[point = 1]))[i]; kw...)
    end
    for g in groups
        if g == :soil_temperature
            soil_panel!("Soil temperature", "°C", :soil_temperature, celsius)
        elseif g == :above_ground
            air_panel!("Air temperature", "°C", :air_temperature, celsius)
            air_panel!("Wind speed", "m/s", :wind_speed, w -> ustrip(u"m/s", w))
            air_panel!("Relative humidity", "%", :relative_humidity, percent)
        elseif g == :radiation
            push!(panels, ("Radiation, $(c.label)", "W/m²", ax -> begin
                series!(ax, out, :global_radiation, r -> ustrip(u"W/m^2", r); label = "solar, horizontal")
                on_slope = slope_solar(out, s)
                on_slope === nothing || lines!(ax, x, on_slope[idx]; label = "solar, on the slope")
                series!(ax, out, :sky_temperature, T -> STEFAN_BOLTZMANN * ustrip(u"K", T)^4; label = "longwave down")
                surface = ustrip.(u"K", collect(out.soil_temperature[point = 1, depth = 1]))[idx]
                lines!(ax, x, SURFACE_EMISSIVITY * STEFAN_BOLTZMANN .* surface .^ 4; label = "longwave up")
            end, identity))
        elseif g == :soil_water
            soil_panel!("Soil moisture", "% by volume", :soil_moisture, percent)
            soil_panel!("Soil humidity", "%", :soil_humidity, percent)
            soil_panel!("Soil water potential, as suction −ψ", "J/kg", :soil_water_potential, suction; yscale = log10)
        elseif g == :snow
            push!(panels, ("Snow depth", "cm", ax -> begin
                for (k, (ck, o)) in enumerate(results)
                    series!(ax, o, :snow_depth, d -> ustrip(u"cm", d); color = Makie.wong_colors()[mod1(k, 7)], label = ck.label)
                end
            end, identity))
        elseif g == :dew_frost
            push!(panels, ("Dew and frost formed at the ground surface, $(c.label)", "kg/m² per hour", ax -> begin
                series!(ax, out, :ground_dew, kg_m2; label = "dew")
                series!(ax, out, :ground_frost, kg_m2; label = "frost")
            end, identity))
        end
    end
    isempty(panels) && return DOM.p("Choose something to plot.")
    fig = Figure(size = (980, 40 + 210 * length(panels)))
    axes = [Axis(fig[k, 1]; title, ylabel, yscale, xlabel = k == length(panels) ? xlabel : "")
            for (k, (title, ylabel, _, yscale)) in enumerate(panels)]
    for (k, (ax, (_, _, draw!, _))) in enumerate(zip(axes, panels))
        draw!(ax)
        labelled_plots = [p for p in ax.scene.plots if haskey(p.attributes, :label) && !isnothing(p.label[])]
        isempty(labelled_plots) || Legend(fig[k, 2], ax; framevisible = false, labelsize = 10, padding = (0, 0, 0, 0))
    end
    # Sunrise and sunset (dashed) on a single day.
    doy = day_shown(v, climatology, s.year)
    if doy !== nothing
        rise, set = sun_hours(s.lat, doy)
        0 < rise < 12 && foreach(ax -> vlines!(ax, [rise, set]; color = (:grey50, 0.7), linestyle = :dash, linewidth = 1), axes)
    end
    if climatology && v.period == :year
        for ax in axes
            ax.xticks = (12:24:287, MONTH_ABBREVIATIONS)
            vlines!(ax, 24:24:263; color = (:black, 0.25), linewidth = 0.5)
        end
    end
    linkxaxes!(axes...)
    return fig
end

# ── The point's terrain ───────────────────────────────────────────────────────────────────────────────────────────

# The terrain a run used: from the SRTM DEM, or CRUCL2's elevation with the sliders' slope, aspect and horizon.
function point_terrain(s::Settings)
    known = lock(() -> get(TERRAIN, (s.lon, s.lat, s.terrain), nothing), CACHE_LOCK)
    known === nothing && return nothing
    return s.terrain ? known : set_terrain!(deepcopy(known), s)
end

function terrain_summary(s::Settings)
    t = point_terrain(s)
    t === nothing && return ""
    h = ustrip.(u"°", collect(first(t.horizon_angles)))
    "Terrain at the point$(s.terrain ? " (from SRTM)" : " (set by hand, elevation from CRUCL2's 10′ grid)"): elevation $(round(ustrip(u"m", first(t.elevation)); digits = 0)) m, " *
        "slope $(round(ustrip(u"°", first(t.slope)); digits = 1))°, aspect $(round(ustrip(u"°", first(t.aspect)); digits = 0))°, " *
        "horizon $(round(minimum(h); digits = 0))–$(round(maximum(h); digits = 0))°"
end

# Terrain maps around the point (elevation, slope, aspect, sky view factor), from the SRTM DEM with the package's own
# terrain code, computed once per point.
const TERRAIN_PANELS = [
    (:elevation, "Elevation"), (:slope, "Slope"), (:aspect, "Aspect"), (:sky_view, "Sky view factor"),
    (:skyline, "Horizon and the sun's path"),
]
const TERRAIN_MAPS = Dict{Tuple{Float64,Float64},Any}()

function terrain_maps(lon, lat; d = 0.02)
    key = (lon, lat)
    known = lock(() -> get(TERRAIN_MAPS, key, nothing), CACHE_LOCK)
    known === nothing || return known
    area = Rasters.Extents.Extent(X = (lon - d, lon + d), Y = (lat - d, lat + d))
    dem = read(crop(Raster(SRTM; extent = area, lazy = true); to = area))
    dem = rebuild(dem; data = [ismissing(z) ? 0.0 : Float64(z) for z in parent(dem)], missingval = nothing)
    grids = MicroclimateMapper.compute_terrain_grids(dem; n_horizon_angles = 32)
    # The outermost cells lack neighbours for slope, aspect and horizons, so they are trimmed.
    inner(a) = a[2:end-1, 2:end-1]
    maps = (; x = collect(lookup(grids.elevation, X))[2:end-1], y = collect(lookup(grids.elevation, Y))[2:end-1],
        elevation = inner(ustrip.(u"m", parent(grids.elevation))), slope = inner(ustrip.(u"°", parent(grids.slope))),
        aspect = inner(ustrip.(u"°", parent(grids.aspect))),
        sky_view = inner(map(MicroclimateMapper._sky_view_from_horizon, parent(grids.horizon_angles))))
    lock(() -> (TERRAIN_MAPS[key] = maps), CACHE_LOCK)
    return maps
end

function terrain_figure(s::Settings, v::View, climatology::Bool, panels)
    t = point_terrain(s)
    t === nothing && return DOM.p("Run first.")
    isempty(panels) && return DOM.p("Choose a terrain panel.")
    fig = Figure(size = (980, 330 * cld(length(panels), 2)))
    maps = any(!=(:skyline), panels) ? (try terrain_maps(s.lon, s.lat) catch e; e end) : nothing
    for (k, p) in enumerate(panels)
        row, col = fldmod1(k, 2)
        slot = fig[row, col]
        if p == :skyline
            # The horizon angle in each direction, the sky below it filled, and the sun's path over it.
            # The axis is centred on the equator-facing direction (south, or north in the southern hemisphere),
            # so the sun's noon crossing is mid-plot rather than split across the ends.
            centre = s.lat < 0 ? 0 : 180
            xpos(az) = mod(az - centre + 180, 360)
            angles = ustrip.(u"°", collect(first(t.horizon_angles)))
            n = length(angles)
            xs = xpos.((0:n-1) .* (360 / n))
            order = sortperm(xs)
            azimuths = vcat(xs[order], 360)
            skyline = vcat(angles[order], angles[order[1]])
            compass = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
            ticks = [compass[mod1(round(Int, (x + centre - 180) / 45) + 1, 8)] for x in 0:45:360]
            ax = Axis(slot[1, 1]; title = "Horizon and the sun's path$(s.terrain ? "" : " (horizon set by hand)")",
                xlabel = "Direction", ylabel = "Elevation angle (°)",
                xticks = (0:45:360, ticks), limits = (0, 360, 0, 90))
            band!(ax, azimuths, zeros(n + 1), skyline; color = (:sienna, 0.5), label = "terrain")
            doy = day_shown(v, climatology, s.year)
            days = doy === nothing ? [(172, "June solstice"), (80, "March equinox"), (355, "December solstice")] :
                [(doy, "the day shown")]
            for (j, (day, label)) in enumerate(days)
                path = sun_path(s.lat, day)
                isempty(path) && continue
                # Break the line where the path crosses the ends of the axis (near the tropics at noon).
                x, y = Float64[], Float64[]
                for (a, e) in path
                    p = xpos(a)
                    !isempty(x) && abs(p - x[end]) > 180 && (push!(x, NaN); push!(y, NaN))
                    push!(x, p); push!(y, e)
                end
                lines!(ax, x, y; color = Makie.wong_colors()[j], label = "sun, $label")
            end
            Legend(slot[1, 2], ax; labelsize = 9, framevisible = false)
        elseif maps isa Exception || maps === nothing
            ax = Axis(slot[1, 1]; title = last(TERRAIN_PANELS[findfirst(q -> first(q) == p, TERRAIN_PANELS)]))
            text!(ax, 0.5, 0.5; text = "DEM unavailable", space = :relative, align = (:center, :center))
        else
            title, colormap, label = p == :elevation ? ("Elevation (SRTM)", :terrain, "m") :
                p == :slope ? ("Slope", :viridis, "°") : p == :aspect ? ("Aspect (clockwise from north)", :phase, "°") :
                ("Sky view factor", :grays, "")
            ax = Axis(slot[1, 1]; title, aspect = DataAspect(), xlabel = "Longitude", ylabel = "Latitude")
            hm = heatmap!(ax, maps.x, maps.y, getproperty(maps, p); colormap,
                colorrange = p == :aspect ? (0, 360) : Makie.automatic)
            Colorbar(slot[1, 2], hm; label)
            scatter!(ax, [s.lon], [s.lat]; color = :white, strokecolor = :black, strokewidth = 1.5, markersize = 12)
        end
    end
    return fig
end

# Every depth and height, every hour, for one data set, written to `path` when the user asks.
function write_csv(path, c::Choice, o)
    soil(layer, f) = [f.(collect(getproperty(o, layer)[point = 1, depth = d])) for d in eachindex(DEPTHS)]
    air(layer, f) = [f.(collect(getproperty(o, layer)[point = 1, height = h])) for h in eachindex(HEIGHTS)]
    header = vcat(["time"],
        ["soil_temperature_$(l)_C" for l in DEPTH_LABELS], ["soil_moisture_$(l)_pct" for l in DEPTH_LABELS],
        ["soil_water_potential_$(l)_J_per_kg" for l in DEPTH_LABELS], ["soil_humidity_$(l)_pct" for l in DEPTH_LABELS],
        ["air_temperature_$(l)_C" for l in HEIGHT_LABELS], ["relative_humidity_$(l)_pct" for l in HEIGHT_LABELS],
        ["snow_depth_cm", "ground_dew_kg_m2", "ground_frost_kg_m2"])
    columns = vcat(soil(:soil_temperature, celsius), soil(:soil_moisture, percent),
        soil(:soil_water_potential, x -> ustrip(u"J/kg", x)), soil(:soil_humidity, percent),
        air(:air_temperature, celsius), air(:relative_humidity, percent),
        [ustrip.(u"cm", collect(o.snow_depth[point = 1])), kg_m2.(collect(o.ground_dew[point = 1])),
         kg_m2.(collect(o.ground_frost[point = 1]))])
    t = collect(lookup(o, Ti))
    open(path, "w") do io
        println(io, join(replace.(header, " " => ""), ","))
        for k in eachindex(t)
            print(io, t[k])
            for col in columns
                print(io, ",", round(col[k]; sigdigits = 5))
            end
            println(io)
        end
    end
    return path
end

file_name(c::Choice, s::Settings) = string(replace(c.code, r"[^A-Za-z0-9]+" => "_"), "_",
    round(s.lon; digits = 3), "_", round(s.lat; digits = 3), c.years === nothing ? "" : "_$(s.year)", ".csv")

# ── The data folder ───────────────────────────────────────────────────────────────────────────────────────────────

default_data_path() = get(ENV, "RASTERDATASOURCES_PATH", joinpath(homedir(), "spatial_data"))

"Point RasterDataSources at `path`, creating it if needed, and describe what is cached there."
function use_data_path!(path::AbstractString)
    path = abspath(expanduser(strip(path)))
    mkpath(path)
    ENV["RASTERDATASOURCES_PATH"] = path
    cached = filter(d -> isdir(joinpath(path, d)), readdir(path))
    return path, isempty(cached) ? "Data folder $path: empty; data sets will be downloaded here when first used." :
        "Data folder $path: already has " * join(cached, ", ") * "."
end

# ── The app ───────────────────────────────────────────────────────────────────────────────────────────────────────

const STYLE = "font-family: sans-serif; font-size: 14px;"
const NOTE = "font-size: 12px; color: #555;"
labelled(text, widget) = DOM.div(DOM.label(text; style = "display:block; font-weight:600; margin-top:8px;"), widget)
checklist(boxes, labels) = DOM.div([DOM.span(b, " ", l, "  "; style = "font-size: 12px; white-space: nowrap;")
    for (b, l) in zip(boxes, labels)]...)

function make_app(data_path::AbstractString = default_data_path())
    App(; title = "MicroclimateMapper.jl") do session
        # Untyped, so that the [lon, lat] array sent from the map's click handler is accepted as it is.
        point = Observable{Any}([-116.545, 33.830])
        mapdiv = DOM.div(; style = "height: 360px; width: 100%; border: 1px solid #ccc;")

        # Inputs, in the order they appear.
        run = Bonito.Button("Run")
        place = Bonito.TextField("Palm Springs, California")
        find = Bonito.Button("Find place")
        data_set = Bonito.Dropdown([c.label for c in CHOICES]; index = 1)
        year = Bonito.NumberInput(2010.0)
        shade = Bonito.Slider(0:5:90; value = 0)
        slope = Bonito.Slider(0:5:90; value = 0)
        aspect = Bonito.Slider(0:15:345; value = 180)
        horizon = Bonito.Slider(0:5:60; value = 0)
        terrain = Bonito.Checkbox(false)
        albedo = Bonito.Slider(0.05:0.05:0.6; value = 0.15)
        soil = Bonito.Dropdown(["example loam", "SoilGrids (texture map)", "SLGA (Australia)"]; index = 1)
        snow = Bonito.Checkbox(false)
        dynamic_moisture = Bonito.Checkbox(false)
        compare = [Bonito.Checkbox(false) for _ in CHOICES]
        data_folder = Bonito.TextField(data_path)
        set_folder = Bonito.Button("Use this folder")
        silo_email = Bonito.TextField(get(ENV, "SILO_EMAIL", ""))
        cds_key = Bonito.TextField(get(ENV, "CDS_API_KEY", ""))
        set_credentials = Bonito.Button("Use these")
        # What to plot.
        group_boxes = [Bonito.Checkbox(g in DEFAULT_GROUPS) for (g, _) in GROUPS]
        terrain_boxes = [Bonito.Checkbox(false) for _ in TERRAIN_PANELS]
        period = Bonito.Dropdown(["the whole year", "a month", "a day"]; index = 1)
        month_choice = Bonito.Dropdown(MONTH_NAMES; index = 1)
        day_choice = Bonito.NumberInput(15.0)
        depth_boxes = [Bonito.Checkbox(l in ("surface", "10 cm", "30 cm", "1 m")) for l in DEPTH_LABELS]
        height_boxes = [Bonito.Checkbox(l in ("1 cm", "30 cm", "2 m")) for l in HEIGHT_LABELS]
        save = Bonito.Button("Save results as CSV")

        status = Observable("Choose a place and a data set, then Run. The first run of a session compiles the model: a minute or two.")
        folder_status = Observable(last(use_data_path!(data_path)))
        credential_status = Observable(haskey(ENV, "CDS_API_KEY") ? "ERA5: CDS key set." : "")
        terrain_text = Observable("")
        figure_slot = Observable{Any}(DOM.div())
        terrain_slot = Observable{Any}(DOM.div())
        script_text = Observable("")
        save_status = Observable("")
        results = Observable{Any}(nothing)
        last_settings = Ref(Settings())

        on(set_folder.value) do _
            try
                folder_status[] = last(use_data_path!(data_folder.value[]))
            catch e
                folder_status[] = "Could not use that folder: " * first(sprint(showerror, e), 200)
            end
        end
        on(set_credentials.value) do _
            isempty(strip(silo_email.value[])) || (ENV["SILO_EMAIL"] = strip(silo_email.value[]))
            era5 = use_cds_era5!(cds_key.value[])
            credential_status[] = join(filter(!isempty, [haskey(ENV, "SILO_EMAIL") ? "SILO email set." : "",
                era5 ? "ERA5 will be read point by point from ECMWF's store (needs `using ZarrDatasets` in this session)." : ""]), " ")
        end
        on(find.value) do _
            try
                g = geocode(place.value[])
                point[] = [g.lon, g.lat]
                # Move the marker to the found place and centre the map on it (a click leaves the view alone). The
                # map's script holds its own copy of `point`, which updates from Julia don't reach, so this is a call.
                lon, lat = g.lon, g.lat
                evaljs(session, js"""(() => {
                    const map = $(mapdiv).leafletMap, marker = $(mapdiv).leafletMarker;
                    marker && marker.setLatLng([$(lat), $(lon)]);
                    map && map.setView([$(lat), $(lon)], Math.max(map.getZoom(), 9));
                })()""")
                status[] = "Found: $(g.display_name)"
            catch e
                status[] = "Place not found: $(sprint(showerror, e))"
            end
        end

        coverage = map(point) do p
            lon, lat = Float64(p[1]), Float64(p[2])
            here = [c.label for c in CHOICES if covers(c, lon, lat)]
            "Lon $(round(lon; digits = 3)), lat $(round(lat; digits = 3)). Data sets covering it: " * join(here, "; ")
        end
        note = map(i -> CHOICES[i].note, data_set.option_index)

        settings() = Settings(; lon = Float64(point[][1]), lat = Float64(point[][2]), year = round(Int, year.value[]),
            shade = shade.value[] / 100, snow = snow.value[], dynamic_moisture = dynamic_moisture.value[], terrain = terrain.value[],
            slope = slope.value[], aspect = aspect.value[], horizon = horizon.value[], albedo = albedo.value[],
            soil = (:example, :SoilGrids, :SLGA)[soil.option_index[]], data_path = ENV["RASTERDATASOURCES_PATH"])
        view() = View((:year, :month, :day)[period.option_index[]], month_choice.option_index[],
            clamp(round(Int, day_choice.value[]), 1, 31))

        shown(boxes) = [k for (k, b) in enumerate(boxes) if b.value[]]
        drawn_once = Ref(false)
        function redraw()
            results[] === nothing && return
            before = status[]
            status[] = drawn_once[] ? "Drawing the figure…" :
                "Drawing the figure… The first time, this compiles the plotting code and takes about a minute."
            sleep(0.1)   # let the message reach the browser before the work starts
            try
                s, v = last_settings[], view()
                figure_slot[] = figure(results[], v, shown(depth_boxes), shown(height_boxes),
                    [GROUPS[k][1] for k in shown(group_boxes)], s)
                panels = [TERRAIN_PANELS[k][1] for k in shown(terrain_boxes)]
                terrain_slot[] = isempty(panels) ? DOM.div() :
                    terrain_figure(s, v, is_climatology(last(first(results[]))), panels)
                drawn_once[] = true
                status[] = before
            catch e
                status[] = "Could not draw the figure: " * first(sprint(showerror, e), 300)
            end
        end

        on(run.value) do _
            s = settings()
            chosen = unique(vcat([CHOICES[data_set.option_index[]]], [CHOICES[k] for k in eachindex(CHOICES) if compare[k].value[]]))
            outside = [c.label for c in chosen if !covers(c, s.lon, s.lat)]
            isempty(outside) || (status[] = "Not covered here: " * join(outside, "; "); return)
            script_text[] = script(chosen, s)
            last_settings[] = s
            Threads.@spawn lock(RUN_LOCK) do
                # The first data set computes the terrain; the rest reuse it and are solved in parallel, one each.
                function attempt(c)
                    try
                        t = @elapsed out = solve_point(c, s)
                        return (c, out, "$(c.label): $(round(t; digits = 1)) s")
                    catch e
                        return (c, nothing, "$(c.label) failed: " * first(sprint(showerror, e), 300))
                    end
                end
                status[] = "Running $(first(chosen).label)… (loading data, then solving)"
                attempts = Any[attempt(first(chosen))]
                if length(chosen) > 1
                    status[] = "Running $(join((c.label for c in chosen[2:end]), ", ")) in parallel…"
                    append!(attempts, fetch.([Threads.@spawn attempt(c) for c in chosen[2:end]]))
                end
                done = Tuple{Choice,Any}[(c, o) for (c, o, _) in attempts if o !== nothing]
                status[] = join(String[m for (_, _, m) in attempts], "  ·  ")
                terrain_text[] = terrain_summary(s)
                isempty(done) || (results[] = done)
            end
        end

        on(_ -> redraw(), results)
        for o in (period.option_index, month_choice.option_index, day_choice.value)
            on(_ -> redraw(), o)
        end
        foreach(b -> on(_ -> redraw(), b.value), vcat(depth_boxes, height_boxes, group_boxes, terrain_boxes))
        on(save.value) do _
            results[] === nothing && (save_status[] = "Run first."; return)
            folder = mkpath(joinpath(ENV["RASTERDATASOURCES_PATH"], "microclimate_outputs"))
            files = [write_csv(joinpath(folder, file_name(c, last_settings[])), c, o) for (c, o) in results[]]
            save_status[] = "Saved: " * join(files, ", ")
        end

        map_js = js"""
        const p = $(point).value;
        const map = L.map($(mapdiv)).setView([p[1], p[0]], 6);
        const osm = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {attribution: '© OpenStreetMap'}).addTo(map);
        const sat = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', {attribution: 'Esri'});
        L.control.layers({"Map": osm, "Satellite": sat}).addTo(map);
        L.control.scale().addTo(map);
        const marker = L.marker([p[1], p[0]]).addTo(map);
        map.on('click', e => { marker.setLatLng(e.latlng); $(point).notify([e.latlng.lng, e.latlng.lat]); });
        // Kept on the element so the Julia-side callbacks below can reach them.
        $(mapdiv).leafletMap = map;
        $(mapdiv).leafletMarker = marker;
        """
        copy = DOM.button("Copy script"; onclick = js"""() => navigator.clipboard.writeText($(script_text).value)""")
        boxes(bs, labels) = DOM.div([DOM.div(b, " ", l; style = "font-size: 12px;") for (b, l) in zip(bs, labels)]...)
        # Soil water and snow rows appear only while dynamic soil moisture or snow is ticked.
        needs = Dict(:soil_water => dynamic_moisture.value, :snow => snow.value)
        row_style(g) = haskey(needs, g) ? map(on -> on ? "font-size: 12px;" : "display: none;", needs[g]) : "font-size: 12px;"
        group_list = DOM.div([DOM.div(b, " ", l; style = row_style(g)) for (b, (g, l)) in zip(group_boxes, GROUPS)]...)

        left = DOM.div(
            DOM.h3("MicroclimateMapper.jl"),
            DOM.div(run; style = "margin-bottom: 8px;"),
            DOM.div(place, find), mapdiv, DOM.p(coverage; style = NOTE),
            labelled("Data set", data_set), DOM.p(note; style = NOTE),
            labelled("Year (for data sets that are not climatologies)", year),
            labelled("Shade (%)", DOM.div(shade, " ", shade.value)),
            labelled("Slope (°)", DOM.div(slope, " ", slope.value)),
            labelled("Aspect (°, clockwise from north: 180 faces south)", DOM.div(aspect, " ", aspect.value)),
            labelled("Horizon angle (°, the same in every direction)", DOM.div(horizon, " ", horizon.value)),
            labelled("Or take elevation, slope, aspect and horizon from the SRTM DEM, about 90 m (slower; the sliders above are then ignored)", terrain),
            labelled("Surface albedo", DOM.div(albedo, " ", albedo.value)),
            labelled("Soil", soil),
            labelled("Snow", snow),
            labelled("Dynamic soil moisture (solved from rainfall; slower)", dynamic_moisture),
            labelled("Also run, to compare", boxes(compare, [c.label for c in CHOICES])),
            labelled("Data folder (RASTERDATASOURCES_PATH)", DOM.div(data_folder, set_folder)),
            DOM.p(folder_status; style = NOTE),
            labelled("SILO email address", silo_email),
            labelled("Copernicus (CDS) API key, for ERA5 and ERA5-Land", cds_key), set_credentials,
            DOM.p(credential_status; style = NOTE),
            DOM.p("Kept in this Julia session only, never written into the script."; style = NOTE);
            style = "width: 380px; min-width: 380px; padding: 10px; overflow: hidden; " * STYLE)
        right = DOM.div(
            DOM.p(status; style = "font-weight: 600;"),
            DOM.p(terrain_text; style = NOTE),
            DOM.div(
                labelled("What to plot", group_list),
                labelled("Terrain around the point", boxes(terrain_boxes, last.(TERRAIN_PANELS)));
                style = "display: flex; gap: 32px;"),
            labelled("Period", DOM.div(period, " ", month_choice, " day ", day_choice)),
            DOM.p("A month or a day of a climatology is its representative day. Sunrise and sunset are marked on a single day.";
                style = NOTE),
            labelled("Soil depths to plot", checklist(depth_boxes, DEPTH_LABELS)),
            labelled("Heights to plot (2 m is the reference height)", checklist(height_boxes, HEIGHT_LABELS)),
            figure_slot, terrain_slot,
            DOM.div(save, " ", DOM.span(save_status; style = NOTE)),
            DOM.h4("Script to reproduce this run"), copy,
            DOM.pre(script_text; style = "font-size: 11px; background: #f6f6f6; padding: 8px; white-space: pre-wrap;");
            style = "flex: 1; min-width: 0; padding: 10px; " * STYLE)
        return DOM.div(LEAFLET_CSS, LEAFLET_JS, DOM.div(left, right; style = "display: flex;"), map_js)
    end
end

function app(; port::Integer = 8080, open::Bool = true, data_path::AbstractString = default_data_path(),
        host::AbstractString = "127.0.0.1", proxy_url = nothing)
    server = proxy_url === nothing ? Bonito.Server(make_app(data_path), host, port) :
        Bonito.Server(make_app(data_path), host, port; proxy_url)
    url = Bonito.online_url(server, "/")
    @info "MicroclimateMapper app running at $url — close the server with `close(server)`."
    open && Bonito.HTTPServer.openurl(url)
    return server
end

end
