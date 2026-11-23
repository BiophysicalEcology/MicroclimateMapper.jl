module MicroclimateMapperZarrExt

# Reads ERA5 from ARCO Zarr stores, for `ARCOTimeSeries` sources.

using Dates
using Zarr
using MicroclimateMapper
using RasterDataSources
using RasterDataSources: ERA5, ECMWF, ARCO
using Rasters
using Rasters: X, Y, Ti, Near
using Rasters.Extents: Extent
using GeoFormatTypes: EPSG

const MM = MicroclimateMapper
const RDS = RasterDataSources

# The store holding `name`: ECMWF has one per group of layers, Google one for all.
_store_source(T::Type{<:ERA5{ECMWF,ARCO}}, name::Symbol) = getraster(T, name)
_store_source(T::Type{<:ERA5}, ::Symbol) = getraster(T)

_headers(::Type{<:ERA5{ECMWF}}) = ["Authorization" => "Bearer $(RDS._cds_token())"]
_headers(::Type) = Pair{String,String}[]

# Opened stores, by url. Chunks are cached on disk in the source's cache folder.
const _STORES = Dict{String,Any}()

function _open(T::Type, name::Symbol)
    source = _store_source(T, name)
    group = get!(_STORES, source.url) do
        remote = Zarr.HTTPStore(source.url; headers = _headers(T))
        zopen(Zarr.CachingStore(remote, Zarr.DirectoryStore(source.cache)); consolidated = true)
    end
    return group, group[layername(T, name)]
end

const _CF_STEPS = Dict("seconds" => Second, "minutes" => Minute, "hours" => Hour, "days" => Day)

# CF time coordinate, e.g. "hours since 1900-01-01 00:00:00": values, step and epoch.
function _time_coordinate(group)
    t = group["time"]
    m = match(r"^(\w+) since (\S+)(?: (\S+))?", t.attrs["units"])
    m === nothing && error("Can't read time units \"$(t.attrs["units"])\"")
    epoch = DateTime(m[3] === nothing ? m[2] : "$(m[2])T$(m[3])")
    return t[:], _CF_STEPS[m[1]], epoch
end

# Index range of the store's time steps from `t0` to `t1`, and their DateTimes.
function _time_indices(group, (t0, t1))
    values, step, epoch = _time_coordinate(group)
    ms = Dates.value(convert(Millisecond, step(1)))
    v0 = Dates.value(DateTime(t0) - epoch) ÷ ms
    v1 = Dates.value(DateTime(t1) - epoch) ÷ ms
    ti = searchsortedfirst(values, v0):searchsortedlast(values, v1)
    isempty(ti) && error("The store has no data from $t0 to $t1")
    return ti, epoch .+ step.(values[ti])
end

function _index_range(coords, (lo, hi), dim)
    i = findall(c -> lo <= c <= hi, coords)
    isempty(i) && error("The store has no $dim coordinates in $lo to $hi")
    return first(i):last(i)
end

function MM._arco_read(T::Type, name::Symbol, area::Extent, times)
    group, array = _open(T, name)
    load_area, restore_dims = MM._native_lon_crop(MM.longitude_convention(T), area)
    lons, lats = group["longitude"][:], group["latitude"][:]
    xi = _index_range(lons, load_area.X, X)
    yi = _index_range(lats, load_area.Y, Y)
    ti, datetimes = _time_indices(group, times)
    data = array[xi, yi, ti]
    spatial = restore_dims((X(lons[xi]), Y(lats[yi])))
    return Raster(data, (spatial..., Ti(datetimes)); crs = EPSG(4326), name)
end

_store_lon(::MM.Longitude180, lon) = lon
_store_lon(::MM.Longitude360, lon) = mod(lon, 360)

# One read per point, of the store's cell nearest to it.
function MM._arco_read_points(T::Type, name::Symbol, points_dim, times)
    group, array = _open(T, name)
    convention = MM.longitude_convention(T)
    lons, lats = group["longitude"][:], group["latitude"][:]
    ti, datetimes = _time_indices(group, times)
    coords = MM.lookup(points_dim)
    out = Matrix{eltype(array)}(undef, length(coords), length(ti))
    for (i, (x, y)) in enumerate(coords)
        xi = argmin(abs.(lons .- _store_lon(convention, x)))
        yi = argmin(abs.(lats .- y))
        out[i, :] = array[xi, yi, ti]
    end
    return Raster(out, (points_dim, Ti(datetimes)); name)
end

end
