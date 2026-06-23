"""
    GeocodeResult

Return type of [`geocode`](@ref). Implements the GeoInterface `PointTrait` so
it can be passed directly as a point in `MicroVectorProblem`. The `extent`
field is always populated from Nominatim's natural bounding box for the place,
expanded by the `buffer` argument (in degrees) supplied to `geocode`.

Fields: `lon`, `lat`, `display_name`, `extent`.
"""
struct GeocodeResult
    lon::Float64
    lat::Float64
    display_name::String
    extent::Extent
end

GeoInterface.isgeometry(::Type{GeocodeResult}) = true
GeoInterface.geomtrait(::GeocodeResult)        = GeoInterface.PointTrait()
GeoInterface.ncoord(::GeoInterface.PointTrait, ::GeocodeResult) = 2
GeoInterface.getcoord(::GeoInterface.PointTrait, r::GeocodeResult, i) =
    i == 1 ? r.lon : r.lat

"""
    load_template(dem_source, site_or_extent) -> Raster

Load a DEM raster over `site_or_extent` for use as the `template` argument of
`MicroRasterProblem`. Accepts a [`GeocodeResult`](@ref) (uses its `extent`
field) or any `Extents.Extent` directly.

# Example
```julia
site     = geocode("Alice Springs, Australia"; buffer = 0.1)
template = load_template(SRTM, site)
prob     = MicroRasterProblem(model; area = site.extent, template, dates, soil_profile)
```
"""
load_template(source, site::GeocodeResult) = load_template(source, site.extent)
function load_template(source, extent::Extent)
    read(crop(Raster(source; extent, lazy = true); to = extent, touches = true))
end

"""
    geocode(place_name; buffer=0.0) -> GeocodeResult

Look up the coordinates of `place_name` using the Nominatim OpenStreetMap
geocoding API. Returns a [`GeocodeResult`](@ref) with fields `lon`, `lat`,
`display_name`, and `extent`.

- `buffer` — degrees added to each side of Nominatim's natural bounding box
  to form `extent`. Default `0.0` uses the place's own bounding box.

`GeocodeResult` implements `GeoInterface.PointTrait`, so it can be passed
directly in the `points` vector of `MicroVectorProblem`. Use `result.extent`
with `MicroRasterProblem`.

# Examples
```julia
# Vector run — single point
site = geocode("Alice Springs, Australia")
prob = MicroVectorProblem(model; points = [site], dates, soil_profile)

# Raster run — natural bbox + 0.1° buffer
site     = geocode("Alice Springs, Australia"; buffer = 0.1)
template = load_template(SRTM, site)
prob     = MicroRasterProblem(model; area = site.extent, template, dates, soil_profile)
```
"""
function geocode(place_name::AbstractString; buffer::Real = 0.0)
    query = URIs.escapeuri(strip(place_name))
    url   = "https://nominatim.openstreetmap.org/search" *
            "?q=$query&format=json&limit=1"

    buf = IOBuffer()
    Downloads.download(url, buf; timeout = 20,
                       headers = ["User-Agent" => "MicroclimateMapper/0.1"])

    seekstart(buf)
    results = JSON3.read(buf)

    isempty(results) && error("No geocoding result for \"$place_name\"")

    r  = first(results)
    bb = r.boundingbox  # [south_lat, north_lat, west_lon, east_lon]

    extent = Extent(
        X = (parse(Float64, bb[3]) - buffer, parse(Float64, bb[4]) + buffer),
        Y = (parse(Float64, bb[1]) - buffer, parse(Float64, bb[2]) + buffer),
    )

    result = GeocodeResult(parse(Float64, r.lon), parse(Float64, r.lat),
                           String(r.display_name), extent)

    @info "Matched: $(result.display_name)"
    return result
end
