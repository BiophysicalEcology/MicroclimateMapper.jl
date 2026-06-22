"""
    geocode(place_name) -> NamedTuple

Look up the coordinates of `place_name` using the Nominatim OpenStreetMap
geocoding API. Returns a `(lon, lat, display_name)` NamedTuple. Pass the
result directly as a point in `MicroVectorProblem` — the `(lon, lat)` fields
are GeoInterface-compatible.

# Example
```julia
point = geocode("Alice Springs, Australia")
prob  = MicroVectorProblem(model; points = [point], dates, soil_profile)
```
"""
function geocode(place_name::AbstractString)
    query = URIs.escapeuri(strip(place_name))
    url   = "https://nominatim.openstreetmap.org/search" *
            "?q=$query&format=json&limit=1"

    buf = IOBuffer()
    Downloads.download(url, buf; timeout = 20,
                       headers = ["User-Agent" => "MicroclimateMapper/0.1"])

    seekstart(buf)
    results = JSON3.read(buf)

    isempty(results) && error("No geocoding result for \"$place_name\"")

    r      = first(results)
    result = (lon          = parse(Float64, r.lon),
              lat          = parse(Float64, r.lat),
              display_name = String(r.display_name))

    @info "Matched: $(result.display_name)"
    return result
end
