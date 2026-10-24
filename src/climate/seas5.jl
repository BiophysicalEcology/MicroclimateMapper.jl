# SEAS5 bindings — ECMWF's seasonal forecast system, via the CDS API
# (RasterDataSources `SEAS5`, needs a CDS API key). A run uses one forecast
# (start month) and one ensemble member, chosen in the source type -- see
# `seas5_forecast`.
#
# Loaded as a daily min/max source like GRIDMET, not sub-daily: CDS only
# serves SEAS5's radiation and precipitation as 24-hourly accumulations, and
# a source has a single native timestep. Daily radiation is what the
# envelope path's clear-sky disaggregation expects anyway. Days are UTC days,
# matching the forecast's 24-hourly lead times.

"""
    SEAS5Forecast{Init,Member}

Weather source for one ECMWF SEAS5 forecast: the forecast started on `Init`
(the 1st of a month) and ensemble member `Member` (0 is the unperturbed
control; 1-50 are perturbed). Build it with [`seas5_forecast`](@ref).
Covers 215 days from `Init`; run `dates` must lie within them.
"""
struct SEAS5Forecast{Init,Member} <: RasterDataSources.RasterDataSource end

"""
    seas5_forecast(init::Date; member = 0) -> Type{<:SEAS5Forecast}

The `weather_source` for the SEAS5 forecast started in `init`'s month, using
ensemble `member` (0-50; 0-24 for the 1981-2016 hindcasts). Run each member
separately for an ensemble of microclimates.

```julia
model = MicroMapModel(; weather_source = seas5_forecast(Date(2026, 9, 1)), ...)
```
"""
function seas5_forecast(init::Date; member::Integer = 0)
    0 <= member <= 50 || throw(ArgumentError("SEAS5 `member` must be in 0:50, got $member"))
    return SEAS5Forecast{Date(year(init), month(init), 1), Int(member)}
end

_seas5_init(::Type{<:SEAS5Forecast{Init}}) where {Init} = Init
_seas5_member(::Type{<:SEAS5Forecast{<:Any,Member}}) where {Member} = Member

const _SEAS5_DAYS = 215
const _SEAS5_DAILY_LEADS = 24:24:(24 * _SEAS5_DAYS)
const _SEAS5_6HOURLY_LEADS = 6:6:(24 * _SEAS5_DAYS)

struct SingleForecast <: Loader end

weather_calendar(::Type{<:SEAS5Forecast}) = Daily()
loader(::Type{<:SEAS5Forecast}) = SingleForecast()
# Half SEAS5's 1° grid step: enough to always take in the nearest grid
# point (they sit at cell centres, e.g. 144.5, -37.5).
points_load_buffer(::Type{<:SEAS5Forecast}) = 0.5

function variables(::Type{<:SEAS5Forecast})
    (
        Variable(Temperature(Maximum()), :mx2t24, u"K"),
        Variable(Temperature(Minimum()), :mn2t24, u"K"),
        Variable(DewpointTemperature(), :d2m, u"K"),
        Variable(WindSpeed(), :ws10, u"m/s"),
        Variable(CloudCover(), :tcc, 1),
        Variable(CloudCover(Minimum()), :tcc_min, 1),
        Variable(CloudCover(Maximum()), :tcc_max, 1),
        Variable(SeaLevelPressure(), :msl, u"Pa"),
        Variable(GlobalRadiation(), :ssrd, u"J/m^2/d"),
        # `:tp` is metres of water per day; × 1000 → kg/m².
        Variable(Rainfall(), :tp, u"kg/m^2", raw -> raw * 1000.0),
        Variable(Elevation(), :orog, u"m"),
    )
end

# Native fields are daily series built from one or more CDS variables:
# `cds` are the variables requested, at `leadtime` hours, and `daily` turns
# their (lon, lat, lead) arrays into one (lon, lat, day) array.
_daily(f, x) = dropdims(f(reshape(x, size(x, 1), size(x, 2), 4, :); dims = 3); dims = 3)
_daily_mean(x) = _daily(mean, x)
_daily_min(x) = _daily(minimum, x)
_daily_max(x) = _daily(maximum, x)
# Accumulated from the forecast start; negative steps are packing noise.
_daily_accumulation(x) = max.(diff(cat(zeros(size(x, 1), size(x, 2), 1), x; dims = 3); dims = 3), 0.0)

const _SEAS5_FIELDS = (
    mx2t24 = (cds = ("maximum_2m_temperature_in_the_last_24_hours",), leadtime = _SEAS5_DAILY_LEADS, daily = identity),
    mn2t24 = (cds = ("minimum_2m_temperature_in_the_last_24_hours",), leadtime = _SEAS5_DAILY_LEADS, daily = identity),
    d2m = (cds = ("2m_dewpoint_temperature",), leadtime = _SEAS5_6HOURLY_LEADS, daily = _daily_mean),
    # Largest 6-hourly speed: the daily maximum the envelope path expects
    ws10 = (cds = ("10m_u_component_of_wind", "10m_v_component_of_wind"), leadtime = _SEAS5_6HOURLY_LEADS,
        daily = (u, v) -> _daily_max(sqrt.(u .^ 2 .+ v .^ 2))),
    tcc = (cds = ("total_cloud_cover",), leadtime = _SEAS5_6HOURLY_LEADS, daily = _daily_mean),
    tcc_min = (cds = ("total_cloud_cover",), leadtime = _SEAS5_6HOURLY_LEADS, daily = _daily_min),
    tcc_max = (cds = ("total_cloud_cover",), leadtime = _SEAS5_6HOURLY_LEADS, daily = _daily_max),
    msl = (cds = ("mean_sea_level_pressure",), leadtime = _SEAS5_6HOURLY_LEADS, daily = _daily_mean),
    ssrd = (cds = ("surface_solar_radiation_downwards",), leadtime = _SEAS5_DAILY_LEADS, daily = _daily_accumulation),
    tp = (cds = ("total_precipitation",), leadtime = _SEAS5_DAILY_LEADS, daily = _daily_accumulation),
)

# Whole-degree box around `area`: CDS returns the grid points (cell centres)
# inside it, so a point-like area still gets its own cell.
function _seas5_request_extent(area::Extent)
    x0, x1 = floor(area.X[1]), ceil(area.X[2])
    y0, y1 = floor(area.Y[1]), ceil(area.Y[2])
    x1 == x0 && (x1 += 1)
    y1 == y0 && (y1 += 1)
    x0, x1 = max(x0, -180.0), min(x1, 180.0)
    y0, y1 = max(y0, -90.0), min(y1, 90.0)
    x1 == x0 && (x0 -= 1)
    y1 == y0 && (y0 -= 1)
    return Extent(X = (x0, x1), Y = (y0, y1))
end

const _SEAS5_LON_DIMS = ("longitude", "lon")
const _SEAS5_LAT_DIMS = ("latitude", "lat")
const _SEAS5_MEMBER_DIMS = ("number", "realization", "member")

# One CDS SEAS5 NetCDF (a single variable) as a (lon, lat, lead) array for
# `member`. Dimensions are matched by name, as CDS writes them
# (longitude, latitude, forecast_period, forecast_reference_time, number).
function _seas5_read(path, member::Int, nlead::Int)
    NCDataset(path) do ds
        gridded = [k for k in keys(ds) if
            any(in(_SEAS5_LON_DIMS), NCDatasets.dimnames(ds[k])) &&
            any(in(_SEAS5_LAT_DIMS), NCDatasets.dimnames(ds[k]))]
        length(gridded) == 1 || error("expected one gridded variable in $path, found $gridded")
        var = ds[only(gridded)]
        dimnames = collect(NCDatasets.dimnames(var))
        lond = dimnames[findfirst(in(_SEAS5_LON_DIMS), dimnames)]
        latd = dimnames[findfirst(in(_SEAS5_LAT_DIMS), dimnames)]
        lons, lats = Float64.(ds[lond][:]), Float64.(ds[latd][:])
        idx = map(dimnames) do d
            d in _SEAS5_MEMBER_DIMS || return Colon()
            m = findfirst(==(member), Int.(ds[d][:]))
            m === nothing && error("SEAS5 ensemble member $member not in $path " *
                "(hindcasts before 2017 have members 0-24)")
            m
        end
        data = var[idx...]
        rest = filter(d -> !(d in _SEAS5_MEMBER_DIMS), dimnames)
        others = [i for (i, d) in enumerate(rest) if d != lond && d != latd]
        data = permutedims(data, [findfirst(==(lond), rest); findfirst(==(latd), rest); others])
        length(data) == length(lons) * length(lats) * nlead || error(
            "SEAS5 file $path has dimensions $(Tuple(rest)) of size $(size(data)), " *
            "expected $nlead lead times")
        values = reshape(Float64[ismissing(x) ? NaN : x for x in data], length(lons), length(lats), nlead)
        return (; lons, lats, values, units = get(var.attrib, "units", ""))
    end
end

# CDS queues each request; a few at a time keeps within its per-user limits.
const _SEAS5_CONCURRENT_REQUESTS = 4

# Daily layers on the full-`years` calendar the other loaders use, NaN outside
# the forecast's 215 days (`_check_weather_dates` keeps runs inside them).
function _load_layers(::SingleForecast, source, fields::Tuple, area::Extent, years)
    init, member = _seas5_init(source), _seas5_member(source)
    extent = _seas5_request_extent(area)
    requests = unique((v, _SEAS5_FIELDS[f].leadtime) for f in fields for v in _SEAS5_FIELDS[f].cds)
    @info "  loading SEAS5 forecast $init, member $member: $(length(requests)) CDS request(s)..."
    paths = Dict(zip(first.(requests), asyncmap(requests; ntasks = _SEAS5_CONCURRENT_REQUESTS) do (v, leadtime)
        getraster(RasterDataSources.SEAS5, v; date = init, extent, leadtime)
    end))
    dates = _daily_date_sequence(years)
    layers = map(fields) do f
        spec = _SEAS5_FIELDS[f]
        reads = map(v -> _seas5_read(paths[v], member, length(spec.leadtime)), spec.cds)
        daily = spec.daily(map(r -> r.values, reads)...)
        (; lons, lats) = first(reads)
        out = fill(NaN, length(lons), length(lats), length(dates))
        for k in axes(daily, 3)
            i = Dates.value(init + Day(k - 1) - first(dates)) + 1
            1 <= i <= length(dates) && (out[:, :, i] .= view(daily, :, :, k))
        end
        Raster(out, (X(lons), Y(lats), Ti(1:length(dates))); crs = EPSG(4326), name = f)
    end
    return NamedTuple{fields}(layers)
end

# Grid elevation from SEAS5's orography, served as geopotential at lead time 0
_load_static(::SingleForecast, _source, ::Tuple{}, _area::Extent) = NamedTuple()
function _load_static(::SingleForecast, source, fields::Tuple, area::Extent)
    extent = _seas5_request_extent(area)
    layers = map(fields) do f
        path = getraster(RasterDataSources.SEAS5, "orography"; date = _seas5_init(source), extent, leadtime = 0)
        (; lons, lats, values, units) = _seas5_read(path, _seas5_member(source), 1)
        geopotential = occursin("s**-2", units) || occursin("s-2", units)
        elevation = geopotential ? values[:, :, 1] ./ ustrip(u"m/s^2", Unitful.gn) : values[:, :, 1]
        Raster(elevation, (X(lons), Y(lats)); crs = EPSG(4326), name = f)
    end
    return NamedTuple{fields}(layers)
end

function _check_weather_dates(source::Type{<:SEAS5Forecast}, dates::AbstractVector{Date})
    init = _seas5_init(source)
    last_day = init + Day(_SEAS5_DAYS - 1)
    (first(dates) >= init && last(dates) <= last_day) || throw(ArgumentError(
        "SEAS5 forecast started $init covers $init to $last_day; run dates " *
        "$(first(dates)) to $(last(dates)) fall outside it"))
    return nothing
end
