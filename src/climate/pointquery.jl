"""
    PointQuery

Loads each layer via `PointDataSources.getpoint` for the point at the centre
of `area`, instead of downloading/cropping a raster grid. Single-point runs
only -- `area` must be point-like (near-zero width) or this throws.
"""
struct PointQuery <: WeatherLoader end

function _load_layers(::PointQuery, source, fields::Tuple, area::Extent, years)
    width_x, width_y = area.X[2] - area.X[1], area.Y[2] - area.Y[1]
    max(width_x, width_y) < 1e-6 || throw(ArgumentError(
        "PointQuery requires a point-like area (got X width $width_x, Y width $width_y) -- " *
        "it only ever loads one location and would silently mis-load a genuine multi-point run"
    ))
    lon, lat = (area.X[1] + area.X[2]) / 2, (area.Y[1] + area.Y[2]) / 2
    date_range = (Date(first(years), 1, 1), Date(last(years), 12, 31))
    layers = map(fields) do name
        @info "  querying $source $name at ($lon, $lat)..."
        nt = PointDataSources.getpoint(source, name; lon, lat, date=date_range)
        Raster(reshape(nt.values, 1, 1, length(nt.values)), (X([lon]), Y([lat]), Ti(nt.times)))
    end
    return NamedTuple{fields}(layers)
end
