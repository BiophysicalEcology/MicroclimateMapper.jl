module FigureHelpers

using CairoMakie
using Markdown
using Unitful
using Rasters
using Rasters: X, Y, Ti, lookup

export figure_axis, markdown_table, celsius, map_panel!, map_figure, mark_sites!, drape_panels, animate_drape,
    series_figure

"""
    figure_axis(xlabel, ylabel; size=(700, 500), kw...)

A `Figure` and `Axis` with minor ticks and grid lines, as used for the figures of the manual.
"""
function figure_axis(xlabel, ylabel; size=(700, 500), kw...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1];
        xlabel, ylabel,
        xminorticksvisible=true, yminorticksvisible=true,
        xminorgridvisible=true, yminorgridvisible=true,
        xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
        kw...,
    )
    return fig, ax
end

"""
    markdown_table(header, rows)

A Markdown table with column names `header` and one row per element of `rows`, with numbers and quantities
rounded to 4 significant digits.
"""
function markdown_table(header, rows)
    lines = ["| " * join(header, " | ") * " |", "|" * repeat(":--|", length(header))]
    append!(lines, ["| " * join(map(_format, row), " | ") * " |" for row in rows])
    return Markdown.parse(join(lines, "\n"))
end

_format(x::AbstractFloat) = string(round(x; sigdigits=4))
_format(x::Quantity) = string(round(unit(x), x; sigdigits=4))
_format(x) = string(x)

"Strip a temperature to a number in °C."
celsius(T) = ustrip(u"°C", T)
celsius(::Missing) = missing

# Cells that are missing or NaN are drawn transparent.
_plain(A) = [ismissing(v) ? NaN : Float64(v isa Quantity ? ustrip(v) : v) for v in A]

"""
    map_panel!(position, raster; title, colormap=:thermal, label="", colorrange=nothing)

A heatmap of an X×Y raster in its own axis at `position` with an equal-aspect geographic frame and a colour bar
beside it. Returns the axis.
"""
function map_panel!(pos, raster; title="", colormap=:thermal, label="", colorrange=nothing, kw...)
    ax = Axis(pos[1, 1]; title, xlabel="Longitude", ylabel="Latitude", aspect=DataAspect(), kw...)
    xs, ys = collect(lookup(raster, X)), collect(lookup(raster, Y))
    values = _plain(parent(raster))
    hm = colorrange === nothing ? heatmap!(ax, xs, ys, values; colormap, nan_color=:transparent) :
         heatmap!(ax, xs, ys, values; colormap, colorrange, nan_color=:transparent)
    Colorbar(pos[1, 2], hm; label)
    return ax
end

"""
    map_figure(panels; ncols=2, size=nothing)

A figure of map panels, each `(raster, title, colormap, label)`.
"""
function map_figure(panels; ncols=2, size=nothing, sites=())
    nrows = cld(length(panels), ncols)
    fig = Figure(; size=something(size, (ncols * 420, nrows * 360)))
    for (k, (raster, title, colormap, label)) in enumerate(panels)
        row, col = fldmod1(k, ncols)
        ax = map_panel!(fig[row, col], raster; title, colormap, label)
        mark_sites!(ax, sites)
    end
    return fig
end

"""
    mark_sites!(ax, sites)

Mark named sites, each `name => (lon, lat)`, with a white dot and a label.
"""
function mark_sites!(ax, sites)
    for (name, (lon, lat)) in sites
        scatter!(ax, [lon], [lat]; color=:white, strokecolor=:black, strokewidth=1, markersize=9)
        text!(ax, lon, lat; text=" " * name, align=(:left, :center), fontsize=11)
    end
    return ax
end

"""
    series_figure(times, series; ylabel, size=(760, 320))

Lines through time, one per `label => values`.
"""
function series_figure(times, series; ylabel="", xlabel="", size=(760, 320), title="")
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; ylabel, xlabel, title)
    for (label, values) in series
        lines!(ax, times, _plain(values); label)
    end
    axislegend(ax; position=:lt, labelsize=10)
    return fig
end

# ── Terrain drapes, from the MicroclimateTalk figures (rafaqz/MicroclimateTalk, src/plot_demos.jl) ──────────

const DRAPE_Z_OFFSET = 5.0

function _drape_clims(raster)
    valid = filter(!isnan, _plain(parent(raster)))
    isempty(valid) && return (0.0, 1.0)
    lo, hi = extrema(valid)
    return lo == hi ? (lo, lo + 1.0) : (lo, hi)
end

function _terrain_base!(ax, xs, ys, zs)
    surface!(ax, xs, ys, zs; color=zs, colormap=:greys, colorrange=extrema(zs), shading=NoShading)
end

"""
    drape_panels(layer, dem, indices, labels; title, colormap=:thermal, label="", ncols=3)

Panels of an X×Y×Ti `layer` at the time `indices`, each draped over the elevation `dem` in 3-D.
"""
function drape_panels(layer, dem, indices, labels; title="", colormap=:thermal, label="", ncols=3,
                      colorrange=_drape_clims(layer))
    xs, ys = collect(lookup(dem, X)), collect(lookup(dem, Y))
    zs = _plain(parent(dem))
    nrows = cld(length(indices), ncols)
    fig = Figure(; size=(ncols * 380, nrows * 320 + 40))
    isempty(title) || Label(fig[0, 1:ncols], title; fontsize=15, halign=:left)
    for (k, i) in enumerate(indices)
        row, col = fldmod1(k, ncols)
        ax = Axis3(fig[row, col]; title=labels[k], xlabel="Lon", ylabel="Lat", zlabel="m",
            azimuth=-π / 4, elevation=π / 8, aspect=(1, 1, 0.35))
        _terrain_base!(ax, xs, ys, zs)
        surface!(ax, xs, ys, zs .+ DRAPE_Z_OFFSET; color=_plain(parent(view(layer; Ti=i))),
            colormap, colorrange, nan_color=:transparent, shading=NoShading)
    end
    Colorbar(fig[1:nrows, ncols + 1]; colormap, colorrange, label)
    return fig
end

"""
    animate_drape(layer, dem, labels, file; title, colormap=Reverse(:ice), label="cm", framerate=4)

Record an X×Y×Ti `layer` draped over `dem`, one frame per time step, to `file`.
"""
function animate_drape(layer, dem, labels, file; title="", colormap=Reverse(:ice), label="cm", framerate=4)
    xs, ys = collect(lookup(dem, X)), collect(lookup(dem, Y))
    zs = _plain(parent(dem))
    colorrange = _drape_clims(layer)
    fig = Figure(; size=(820, 720))
    title_obs = Observable("$title · $(labels[1])")
    ax = Axis3(fig[1, 1]; title=title_obs, xlabel="Lon", ylabel="Lat", zlabel="m",
        azimuth=-π / 4, elevation=π / 8, aspect=(1, 1, 0.35))
    _terrain_base!(ax, xs, ys, zs)
    color = Observable(_plain(parent(view(layer; Ti=1))))
    surface!(ax, xs, ys, zs .+ DRAPE_Z_OFFSET; color, colormap, colorrange, nan_color=:transparent,
        shading=NoShading)
    Colorbar(fig[1, 2]; colormap, colorrange, label)
    record(fig, file, 1:size(layer, Ti); framerate) do k
        color[] = _plain(parent(view(layer; Ti=k)))
        title_obs[] = "$title · $(labels[k])"
    end
    return file
end

end
