"""
    geocode(place_name) -> (lon, lat)

Look up the longitude and latitude of `place_name` using the Nominatim
OpenStreetMap geocoding API. Returns a `(lon, lat)` tuple suitable for use
as a point in `MicroVectorProblem`.

# Example
```julia
point = geocode("Alice Springs, Australia")
prob  = MicroVectorProblem(model; points = [point], dates, soil_profile)
```
"""
function geocode(place_name::AbstractString)
    encoded = replace(strip(place_name), r"\s+" => "%20")
    url = "https://nominatim.openstreetmap.org/search?q=$(encoded)&format=json&limit=1"
    buf = IOBuffer()
    Downloads.download(url, buf; headers = ["User-Agent" => "Julia/MicroclimateMapper"])
    results = JSON3.read(String(take!(buf)))
    isempty(results) && error("geocode: no results for \"$place_name\"")
    r = first(results)
    lon, lat = parse(Float64, r.lon), parse(Float64, r.lat)
    @info "geocode: $(r.display_name)"
    return (lon, lat)
end
