using CairoMakie
using Rasters
using Rasters: X, Y, Ti, lookup
using Unitful
using Dates

const DRAPE_Z_OFFSET = 5.0

_skip_cell(v, mv) =
    ismissing(v) || (!ismissing(mv) && v == mv) || (v isa AbstractFloat && isnan(v))

function drape_color(slice)
    mv = missingval(slice)
    A  = parent(slice)
    return [_skip_cell(v, mv) ? NaN : Float64(v) for v in A]
end

function drape_clims(raster)
    A  = parent(raster)
    mv = missingval(raster)
    valid = Float64[]
    for v in A
        _skip_cell(v, mv) && continue
        push!(valid, Float64(v))
    end
    isempty(valid) && return nothing
    lo, hi = minimum(valid), maximum(valid)
    return lo == hi ? (lo, lo + 1.0) : (lo, hi)
end

function _add_terrain_base!(ax, xs, ys, zs)
    surface!(ax, xs, ys, zs;
        color = zs, colormap = :greys,
        colorrange = (minimum(zs), maximum(zs)),
        shading = NoShading)
end

# Six-panel drape of an X-Y-Ti layer onto a DEM (e.g. soil temperature at 5 cm
# at six hours of one day). `panel_hours` picks the hours-of-day; matching
# Ti indices are resolved against the first calendar day in `layer`.
function plot_raster_temperature_panels(layer, dem;
        panel_hours    = (0, 6, 9, 12, 15, 18),
        panel_labels   = ("Midnight", "Dawn", "Mid-morning",
                          "Midday", "Mid-afternoon", "Dusk"),
        variable_label = "Soil temperature at 5 cm (°C)",
        unit_label     = "°C",
        cmap           = Reverse(:RdYlBu),
        ncols          = 3,
    )
    ti = lookup(layer, Ti)
    day_one = Dates.Date(first(ti))
    panel_indices = [findfirst(t -> Dates.Date(t) == day_one && Dates.hour(t) == h, ti)
                     for h in panel_hours]
    @assert !any(isnothing, panel_indices) "Could not resolve all panel hours in Ti lookup"
    clims = drape_clims(layer)
    clims === nothing && error("layer has no valid cells")
    xs = collect(lookup(dem, X))
    ys = collect(lookup(dem, Y))
    zs = Float64.(parent(dem))
    zs_drape = zs .+ DRAPE_Z_OFFSET
    nrows = cld(length(panel_indices), ncols)
    fig = Figure(size = (ncols * 460, nrows * 380 + 80))
    Label(fig[0, 1:ncols], variable_label;
        fontsize = 15, halign = :left, padding = (10, 0, 0, 0))
    for k in eachindex(panel_indices)
        row, col = fldmod1(k, ncols)
        ax = Axis3(fig[row, col];
            title    = panel_labels[k],
            xlabel   = "Lon", ylabel = "Lat", zlabel = "m",
            azimuth  = -π/4, elevation = π/8,
            aspect   = (1, 1, 0.35),
        )
        _add_terrain_base!(ax, xs, ys, zs)
        surface!(ax, xs, ys, zs_drape;
            color      = drape_color(view(layer; Ti = panel_indices[k])),
            colormap   = cmap, colorrange = clims,
            nan_color  = :transparent, shading = NoShading,
        )
    end
    Colorbar(fig[1:nrows, ncols + 1]; colormap = cmap, colorrange = clims, label = unit_label)
    return fig
end

# Snow-depth line plot at noon over the run window. One line per point.
function plot_vector_snow_depth(output, point_labels)
    snow_depth = output.snow_depth
    ti         = lookup(snow_depth, Ti)
    noon_idx   = findall(t -> Dates.hour(t) == 12, ti)
    isempty(noon_idx) && (noon_idx = collect(eachindex(ti)))
    noon_dates = collect(ti[noon_idx])
    colors     = [RGBf(0.30, 0.30, 0.30),
                  RGBf(0.55, 0.55, 0.55),
                  RGBf(0.63, 0.00, 0.00)]
    fig  = Figure(size = (1000, 500))
    yr   = Dates.year(first(ti))
    mfst = Dates.monthname(Dates.month(first(ti)))
    mlst = Dates.monthname(Dates.month(last(ti)))
    title = mfst == mlst ? "Snow depth — $(mfst) $(yr)" :
                           "Snow depth — $(mfst)–$(mlst) $(yr)"
    axis = Axis(fig[1, 1]; title, xlabel = "Date", ylabel = "Snow depth (cm)")
    for (i, name) in enumerate(point_labels)
        depths_cm = ustrip.(u"cm",
            parent(view(snow_depth; point = i, Ti = noon_idx)))
        lines!(axis, noon_dates, depths_cm;
            label = name, color = colors[mod1(i, length(colors))], linewidth = 1.4)
    end
    axislegend(axis; position = :rt, framevisible = false)
    return fig
end
