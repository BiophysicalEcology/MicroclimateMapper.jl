# Generic per-pixel weather assembler.
#
# A source declares `weather_loader(source)` (how its files are organised)
# and `weather_variables(source)` (which canonical variables its layers
# provide, with native unit/transform). `_load_weather`/`assemble_weather!`
# then load, merge in any fallback layers, and run derivations for whatever
# canonical variables aren't provided natively.
#
# New sources just need those two methods, plus fallback_source/
# fallback_layers/primary_layers if they need a baseline for missing fields.

# ---------------------------------------------------------------------------
# Variable map: canonical → native (per source)
# ---------------------------------------------------------------------------

"""
    WeatherVariable{Name, Field, U, T}(unit, transform)

Declares that raster layer `Field` provides canonical variable `Name`,
with native `unit` (auto-converted on assignment) and optional `transform`
for non-unit scaling. `Name`/`Field` are type parameters, so a tuple of
these gives the compiler one specialised path per entry.
"""
struct WeatherVariable{Name, Field, U, T}
    unit::U
    transform::T
end
WeatherVariable(name::Symbol, field::Symbol, unit = 1, transform = identity) =
    WeatherVariable{name, field, typeof(unit), typeof(transform)}(unit, transform)

@inline canonical_name(::WeatherVariable{Name}) where {Name} = Name
@inline native_field(::WeatherVariable{<:Any, Field}) where {Field} = Field

"""
    weather_variables(::Type{<:RasterDataSource}) -> Tuple{WeatherVariable, …}

Tuple of `WeatherVariable`s declaring which canonical physical variables
the source natively provides. The set of derivations that run is
`weather_derivations(source)` minus those whose output appears in this
native set.
"""
function weather_variables end

# ---------------------------------------------------------------------------
# Temporal-resolution trait
# ---------------------------------------------------------------------------

"""
    TemporalResolution

Trait for a source's time step. Drives per-pixel buffer sizing and which
env-minmax struct gets built (`MonthlyMinMaxEnvironment` vs
`DailyMinMaxEnvironment`).
"""
abstract type TemporalResolution end

"""
    MonthlyResolution

12 steps/year. TerraClimate, CHELSA, WorldClim. Uses `MonthlyMinMaxEnvironment`
and `Microclimate.DEFAULT_DAYS` (mid-month days) for solar geometry.
"""
struct MonthlyResolution <: TemporalResolution end

"""
    DailyResolution

365/366 steps/year (leap-aware). NCEP, AWAP, ERA5 daily aggregates. Uses
`DailyMinMaxEnvironment`; solver runs consecutive-day mode (soil state
carries day-to-day).
"""
struct DailyResolution <: TemporalResolution end

"""
    HourlyResolution

8760/8784 steps/year. ERA5 and similar. No min/max envelope
(`environment_minmax = nothing`) — hourly values feed `environment_hourly`
directly; `environment_daily` is filled by aggregating them.
"""
struct HourlyResolution <: TemporalResolution end

"""
    SixHourlyResolution

Native 6-hourly files (4/day), interpolated to hourly before the solver
sees them. Used by NCEP. Like `HourlyResolution`: no min/max envelope,
consecutive-day mode. Future sub-daily sources should follow the same
pattern with their own resolution trait.
"""
struct SixHourlyResolution <: TemporalResolution end

"""
    temporal_resolution(::Type{<:RasterDataSource}) -> TemporalResolution

The time step a source's `_load_weather` produces. Default
`MonthlyResolution()`; daily/hourly sources must override.
"""
@inline temporal_resolution(::Type) = MonthlyResolution()

# Number of steps/year in the *native* files — drives staging-buffer sizes.
# For sources that don't distinguish native from output resolution, these equal steps_per_year.
@inline _native_steps_per_year(::MonthlyResolution) = 12
@inline _native_steps_per_year(::DailyResolution) = 365
@inline _native_steps_per_year(::HourlyResolution) = 8760
@inline _native_steps_per_year(::SixHourlyResolution) = 1460   # 4 × 365

@inline _minmax_env_type(::MonthlyResolution) = MonthlyMinMaxEnvironment
@inline _minmax_env_type(::DailyResolution) = DailyMinMaxEnvironment

# `years` is the actual calendar-year span (e.g. `2000:2004`), not just a
# count — daily/hourly/sub-daily resolutions need the real per-year day
# count (365 or 366) to produce a correctly-sized day-of-year vector.
@inline _days_of_year(::MonthlyResolution, years) = repeat(Microclimate.DEFAULT_DAYS, length(years))
@inline _days_of_year(::Union{DailyResolution, HourlyResolution, SixHourlyResolution}, years) =
    # `init=Int[]` forces a real `vcat` even for single-year `years` —
    # otherwise `reduce` short-circuits to a bare unconverted UnitRange,
    # which fails SolarRadiation.solar_radiation!'s `Vector{<:Real}` check.
    reduce(vcat, (1:Dates.daysinyear(y) for y in years); init = Int[])

# Per-year step counts (365 or 366, × the sub-daily factor), used to build
# `year_offsets` once per worker (see `_allocate_weather_buffers`).
@inline _year_lengths(::MonthlyResolution, years) = fill(12, length(years))
@inline _year_lengths(::DailyResolution, years) = Dates.daysinyear.(years)
@inline _year_lengths(::HourlyResolution, years) = Dates.daysinyear.(years)
@inline _year_lengths(::SixHourlyResolution, years) = Dates.daysinyear.(years)

# ---------------------------------------------------------------------------
# Date-range helpers
# ---------------------------------------------------------------------------

# Convert an old-style integer years range to a full-calendar Date range.
_years_to_dates(years::AbstractRange{Int}) =
    Date(first(years), 1, 1):Day(1):Date(last(years), 12, 31)

# Derive the integer year span that must be loaded from weather sources to
# cover all dates.
_years_from_dates(d::Date)                 = year(d):year(d)
_years_from_dates(r::AbstractRange{Date})  = minimum(year, r):maximum(year, r)

# Cumulative day offset at which year `years_v[yi]` begins. Setup-time
# only, so allocation doesn't matter here.
function _cumulative_year_days(years_v::Vector{Int}, yi::Int)
    total = 0
    @inbounds for i in 1:(yi - 1)
        total += Dates.daysinyear(years_v[i])
    end
    return total
end

# Normalise user-supplied dates to a sorted Vector{Date}. Feb 29 is kept —
# leap days are simulated like any other day.
_normalise_dates(d::Date) = [d]
function _normalise_dates(r::AbstractRange{Date})
    v = sort(collect(r))
    isempty(v) && error("No simulation dates supplied.")
    return v
end

# 1-based positional Ti range into the full-years weather stack, plus the
# solver's day-of-year vector. `dates_vec` must already be sorted.
# Daily/hourly/sub-daily variants use real per-year day counts, so leap
# years correctly contribute 366 days.
function _ti_range_for_dates(::MonthlyResolution, years, dates_vec::Vector{Date})
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = (yi_start - 1) * 12 + month(start_d)
    ti_end   = (yi_end   - 1) * 12 + month(end_d)
    days_doy = Int[]
    d = Date(year(start_d), month(start_d), 1)
    stop = Date(year(end_d), month(end_d), 1)
    while d <= stop
        push!(days_doy, Microclimate.DEFAULT_DAYS[month(d)])
        d += Month(1)
    end
    return ti_start, ti_end, days_doy
end

function _ti_range_for_dates(::DailyResolution, years, dates_vec::Vector{Date})
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = _cumulative_year_days(years_v, yi_start) + dayofyear(start_d)
    ti_end   = _cumulative_year_days(years_v, yi_end)   + dayofyear(end_d)
    return ti_start, ti_end, dayofyear.(dates_vec)
end

function _ti_range_for_dates(::HourlyResolution, years, dates_vec::Vector{Date})
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = _cumulative_year_days(years_v, yi_start) * 24 + (dayofyear(start_d) - 1) * 24 + 1
    ti_end   = _cumulative_year_days(years_v, yi_end)   * 24 + dayofyear(end_d) * 24
    return ti_start, ti_end, dayofyear.(dates_vec)
end

# Native stack has 4 Ti steps/day, so the Ti slice is in 6h units.
# The returned day-of-year vector has one entry per calendar day because
# solar geometry and the derivation chain operate at daily granularity
# regardless of sub-daily native resolution.
function _ti_range_for_dates(::SixHourlyResolution, years, dates_vec::Vector{Date})
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = _cumulative_year_days(years_v, yi_start) * 4 + (dayofyear(start_d) - 1) * 4 + 1
    ti_end   = _cumulative_year_days(years_v, yi_end)   * 4 + dayofyear(end_d) * 4
    return ti_start, ti_end, dayofyear.(dates_vec)
end

# One "anchor" Date per solver step — used to carry the year when building
# the output DateTime axis. Monthly: 1st of each selected month. Daily/hourly:
# the dates_vec itself (already one per calendar day).
function _step_anchor_dates(::MonthlyResolution, dates_vec::Vector{Date})
    seen = Set{Tuple{Int,Int}}()
    result = Date[]
    for d in sort(dates_vec)
        k = (year(d), month(d))
        k in seen && continue
        push!(seen, k)
        push!(result, Date(year(d), month(d), 1))
    end
    return result
end
_step_anchor_dates(::TemporalResolution, dates_vec::Vector{Date}) = dates_vec

# Monthly forcing → independent representative days (Fortran "monthly mode");
# daily/hourly/sub-daily forcing → continuous run with state carrying day-to-day.
@inline _time_mode(::MonthlyResolution) = Microclimate.NonConsecutiveDayMode()
@inline _time_mode(::DailyResolution) = Microclimate.ConsecutiveDayMode()
@inline _time_mode(::HourlyResolution) = Microclimate.ConsecutiveDayMode()
@inline _time_mode(::SixHourlyResolution) = Microclimate.ConsecutiveDayMode()

# ---------------------------------------------------------------------------
# Loader traits
# ---------------------------------------------------------------------------

"""
    WeatherLoader

Trait for how a source's raster files are organised.
`_load_layers(::WeatherLoader, source, fields, area, years)` turns
source-field symbols into a `NamedTuple` of 3-D rasters `(X, Y, Ti)`.
"""
abstract type WeatherLoader end

"""
    YearlyTimeSeries

One file per year, all 12 months along `Ti`. TerraClimate.
`getraster(source, name; date = Date(year))`.
"""
struct YearlyTimeSeries <: WeatherLoader end

"""
    MonthlyTimeSeries

One file per real calendar month per layer. BARRA.
`getraster(source, name; date = Date(year, month))`.
"""
struct MonthlyTimeSeries <: WeatherLoader end

"""
    MonthlyClimatology

12 files per layer, one per month, fixed climatology (tiled across
`years`). CHELSA{Climate}, WorldClim{Climate}.
`getraster(source, name; month)`.
"""
struct MonthlyClimatology <: WeatherLoader end

"""
    FutureMonthlyClimatology

Like `MonthlyClimatology`, with an extra `date` selecting the projected
window. CHELSA{Future{Climate}}.
`getraster(source, name; date = future_date, month)`.
"""
struct FutureMonthlyClimatology <: WeatherLoader end

"""
    MultiBandFutureClimatology

One multi-band file per layer (12 months as GeoTIFF bands) for a future
period. WorldClim{Future{Climate}}.
`getraster(source, name; date = future_date)`.
"""
struct MultiBandFutureClimatology <: WeatherLoader end

"""
    SingleFileBands

All variables as named bands in one file. `getraster(source)` returns the
path; `fields` are read together as a `RasterStack` and the native time
dimension is tiled `nyears` times. CRUCL2 (12 monthly steps).
"""
struct SingleFileBands <: WeatherLoader end

"""
    DailyFiles

One 2-D raster file per calendar day per layer. AWAP and similar.
`getraster(source, name; date = day)`.
"""
struct DailyFiles <: WeatherLoader end

"""
    HourlyZarrStore

Single remote Zarr store, every layer a 3-D `(time, lat, lon)` variable.
ERA5 via ARCO-ERA5. `getraster(source)` returns a `CachedCloudSource`
whose `url` is opened once; layers are read by their long name
(`RasterDataSources.layername(source, sym)`).
"""
struct HourlyZarrStore <: WeatherLoader end

"""
    weather_loader(::Type{<:RasterDataSource}) -> WeatherLoader

How this source's files are organised.
"""
function weather_loader end

# ---------------------------------------------------------------------------
# Per-source extension points (with defaults)
# ---------------------------------------------------------------------------

"""
    primary_layers(source) -> Tuple{Symbol, …}

Source-native raster layer names loaded via `weather_loader(source)`.
Defaults to the unique `native_field`s from `weather_variables(source)`.
Override when the primary file set is a strict subset of the declared
variables (the rest coming from a baseline).
"""
@inline primary_layers(source) = _unique_native_fields(weather_variables(source))

"""
    fallback_source(source) -> source or Nothing

Baseline source to fall back on for variables `source` can't provide
itself (e.g. CHELSA{Climate} for CHELSA Future). `nothing` = no fallback.
"""
@inline fallback_source(::Type) = nothing

"""
    fallback_layers(source) -> Tuple{Symbol, …}

Source-field names to load from `fallback_source` and merge in.
"""
@inline fallback_layers(::Type) = ()

"""
    _extra_getraster_kwargs(source) -> NamedTuple

Extra kwargs splatted into every `getraster(source, name; …)` call.
"""
@inline _extra_getraster_kwargs(::Type) = (;)

"""
    weather_grid_elevation(source, weather, I) -> Quantity or nothing

Weather grid's elevation at spatial index `I` (e.g. `450.0u"m"`), or
`nothing` if the source carries none — skips lapse correction
(`grid_elevation = site.elevation`, Δz = 0). Sources with their own
elevation layer (e.g. CRUCL2's `elv`) should override this.
"""
@inline weather_grid_elevation(::Type, _weather, _I) = nothing

# ---------------------------------------------------------------------------
# Generic _load_weather
# ---------------------------------------------------------------------------

"""
    _post_load_stack!(source, stack, years) -> RasterStack

Post-processing hook run after every layer of `stack` is loaded. Default
no-op. `years` is the real calendar-year span, so overrides can compute a
leap-year-aware expected length via `Dates.daysinyear`.

Override when a source's files carry sub-daily Ti the declared
`temporal_resolution` doesn't expect — e.g. `NCEP{SurfaceFlux}` in daily
mode stores 6-hourly data needing averaging to daily.
"""
_post_load_stack!(::Type, stack, _years) = stack

# Average a sub-daily Ti layer down by `factor` steps.
function _aggregate_ti_to_daily(layer::AbstractRaster, factor::Int)
    return Rasters.aggregate(mean, layer, (Ti(factor),))
end

"""
    _load_weather(source, area, years) -> RasterStack

Loads every layer the source contributes to the canonical variable map.
Sources with a `fallback_source` get primary layers from themselves and
fallback layers from the baseline, merged into one `RasterStack`.
Dispatch is via `weather_loader(source)`; `_post_load_stack!` runs after.
"""
function _load_weather(source::Type, area::Extent, years)
    primary_stack = _load_layers(weather_loader(source), source,
                                 primary_layers(source), area, years)
    baseline = fallback_source(source)
    stack = if baseline === nothing
        RasterStack(primary_stack)
    else
        fallback_stack = _load_layers(weather_loader(baseline), baseline,
                                      fallback_layers(source), area, years)
        fallback_stack = _match_fallback_resolution(source, baseline, fallback_stack,
                                                     primary_stack, years)
        RasterStack(merge(primary_stack, fallback_stack))
    end
    # Source files declare a `missingval`, so loaded eltypes are
    # `Union{Missing, T}`. The per-pixel reader does `value * unit`,
    # which throws `convert(Missing, Quantity)` on any masked cell.
    # `replace_missing(..., NaN)` strips the Missing union and lets
    # NaN propagate visibly if any cell is genuinely masked.
    stack = Rasters.replace_missing(stack, NaN)
    return _post_load_stack!(source, stack, years)
end

# Align a fallback source's layers onto `source`'s grid/resolution when they
# differ — e.g. CRUCL2's monthly wind climatology backfilling daily AWAP/SILO.
# `template` supplies the target X/Y grid and Ti axis. Same resolution is a
# no-op; unhandled combinations raise clearly.
_match_fallback_resolution(source::Type, baseline::Type, layers::NamedTuple, template::NamedTuple, years) =
    _match_fallback_resolution(temporal_resolution(source), temporal_resolution(baseline), layers, template, years)
_match_fallback_resolution(::T, ::T, layers::NamedTuple, _template, _years) where {T <: TemporalResolution} = layers
function _match_fallback_resolution(::DailyResolution, ::MonthlyResolution, layers::NamedTuple,
                                    template::NamedTuple, years)
    ref = first(values(template))
    names = keys(layers)
    return NamedTuple{names}(map(l -> _monthly_to_daily(l, ref, years), values(layers)))
end
_match_fallback_resolution(target::TemporalResolution, base::TemporalResolution, ::NamedTuple, _template, _years) =
    error("No resolution-matching path from $(typeof(base)) to $(typeof(target)) for a weather fallback")

# Resample a monthly climatology onto `ref`'s grid, then repeat each of its
# 12*nyears slices across every real day of that month (leap-aware), reusing
# `ref`'s own Ti axis. Setup-time only.
function _monthly_to_daily(layer::AbstractRaster, ref::AbstractRaster, years)
    resampled = Rasters.resample(layer; to = dims(ref, (X, Y)))
    data = parent(resampled)
    ti = dims(ref, Ti)
    out = similar(data, size(data, 1), size(data, 2), length(ti))
    years_v = collect(years)
    d = 0
    for (yi, y) in enumerate(years_v), m in 1:12
        @views month_slice = data[:, :, (yi - 1) * 12 + m]
        for _ in 1:Dates.daysinmonth(y, m)
            d += 1
            out[:, :, d] .= month_slice
        end
    end
    return Raster(out, (dims(resampled)[1:2]..., ti); crs = crs(resampled))
end

# Resample a daily layer onto `ref`'s grid, then place each day's value at
# its first hour (midnight), zero elsewhere — a single daily event on `ref`'s
# hourly Ti axis. Setup-time only.
function _daily_to_hourly_midnight(layer::AbstractRaster, ref::AbstractRaster)
    resampled = Rasters.resample(layer; to = dims(ref, (X, Y)))
    data = parent(resampled)
    ti = dims(ref, Ti)
    out = zeros(eltype(data), size(data, 1), size(data, 2), length(ti))
    for d in 1:size(data, 3)
        out[:, :, (d - 1) * 24 + 1] .= @view data[:, :, d]
    end
    return Raster(out, (dims(resampled)[1:2]..., ti); crs = crs(resampled))
end

# Resample a static layer onto `ref`'s grid and tile across `ref`'s Ti axis
# (CRUCL2 `:elv`, GRIDMET `:elev` convention). Setup-time only.
function _static_to_ti(layer::AbstractRaster, ref::AbstractRaster)
    resampled = Rasters.resample(layer; to = dims(ref, (X, Y)))
    data = parent(resampled)
    ti = dims(ref, Ti)
    tiled = repeat(reshape(data, size(data)..., 1); outer = (1, 1, length(ti)))
    return Raster(tiled, (dims(resampled)[1:2]..., ti); crs = crs(resampled))
end

# ---------------------------------------------------------------------------
# Per-loader file-reading
# ---------------------------------------------------------------------------

# Fetches `fetch(name, item)` for every (name, item) pair over `fields ×
# items` on worker threads; each opens its own file handle, so this is safe
# even though NetCDF/GDAL reads aren't internally threaded (~3x speedup on
# BARRA's network archive, values identical to sequential). A one-off HDF5
# "wrong B-tree signature" corruption during a full-scale run traced to an
# external process overwriting those files concurrently, not a thread bug here.
function _parallel_fetch(fetch, fields::Tuple, items)
    nf, ni = length(fields), length(items)
    rasters = Matrix{Any}(undef, nf, ni)
    Threads.@threads for idx in 1:(nf * ni)
        fi, ii = fldmod1(idx, ni)
        rasters[fi, ii] = fetch(fields[fi], items[ii])
    end
    return rasters
end

function _load_layers(::YearlyTimeSeries, source, fields::Tuple, area::Extent, years)
    extras = _extra_getraster_kwargs(source)
    rasters = _parallel_fetch(fields, years) do name, yr
        @info "  loading $source $name $yr..."
        path = getraster(source, name; date = Date(yr), extras...)
        read(crop(Raster(path; lazy = true); to = area, touches = true))
    end
    layers = map(fi -> cat(rasters[fi, :]...; dims = Ti), 1:length(fields))
    return NamedTuple{fields}(Tuple(layers))
end

function _load_layers(::MonthlyTimeSeries, source, fields::Tuple, area::Extent, years)
    extras = _extra_getraster_kwargs(source)
    year_months = [Date(y, m) for y in years for m in 1:12]
    rasters = _parallel_fetch(fields, year_months) do name, d
        @info "  loading $source $name $(year(d))-$(month(d))..."
        path = getraster(source, name; date = d, extras...)
        read(crop(Raster(path; name, lazy = true); to = area, touches = true))
    end
    layers = map(fi -> cat(rasters[fi, :]...; dims = Ti), 1:length(fields))
    return NamedTuple{fields}(Tuple(layers))
end

function _load_layers(::MonthlyClimatology, source, fields::Tuple, area::Extent, years)
    nyears = length(years)
    layers = map(fields) do name
        @info "  loading $source $name (monthly climatology)..."
        _build_monthly_climatology(area, nyears) do month
            getraster(source, name; month)
        end
    end
    return NamedTuple{fields}(layers)
end

function _load_layers(::FutureMonthlyClimatology, source, fields::Tuple, area::Extent, years)
    nyears = length(years)
    future_date = Date(first(years))
    layers = map(fields) do name
        @info "  loading $source $name (future monthly climatology, $future_date)..."
        _build_monthly_climatology(area, nyears) do month
            getraster(source, name; date = future_date, month)
        end
    end
    return NamedTuple{fields}(layers)
end

function _load_layers(::MultiBandFutureClimatology, source, fields::Tuple,
                      area::Extent, years)
    nyears = length(years)
    future_date = Date(first(years))
    layers = map(fields) do name
        @info "  loading $source $name ($future_date)..."
        path = getraster(source, name; date = future_date)
        full = read(crop(Raster(path; lazy = true); to = area, touches = true))
        # File holds 12 months along the 3rd dim (Band or Ti, depending on
        # the format). Pull out the data array, tile across years, and re-
        # wrap with a proper `Ti` axis.
        data = parent(full)
        size(data, 3) == 12 || error(
            "$source: expected a 12-month multi-band file for $name, got $(size(data, 3)) bands")
        tiled = nyears == 1 ? data : repeat(data; outer = (1, 1, nyears))
        spatial_dims = dims(full)[1:2]
        Raster(tiled, (spatial_dims..., Ti(1:(12 * nyears))); crs = crs(full))
    end
    return NamedTuple{fields}(layers)
end

function _load_layers(::SingleFileBands, source, fields::Tuple, area::Extent, years)
    nyears = length(years)
    @info "  loading $source $(fields)..."
    path = getraster(source)
    raw = read(crop(RasterStack(path; name = fields, lazy = true); to = area, touches = true))
    # Infer native step count from the first 3-D layer (some layers, e.g. CRUCL2 :elv, are static 2-D).
    nsteps = 1
    for name in fields
        ndims(raw[name]) >= 3 && (nsteps = size(parent(raw[name]), 3); break)
    end
    ntotal = nsteps * nyears
    layers = map(fields) do name
        lyr  = raw[name]
        data = parent(lyr)
        tiled = if ndims(data) == 2
            repeat(reshape(data, size(data)..., 1); outer = (1, 1, ntotal))
        else
            nyears == 1 ? data : repeat(data; outer = (1, 1, nyears))
        end
        Raster(tiled, (dims(lyr)[1:2]..., Ti(1:ntotal)); crs = crs(lyr))
    end
    return NamedTuple{fields}(layers)
end

function _load_layers(::DailyFiles, source, fields::Tuple, area::Extent, years)
    extras = _extra_getraster_kwargs(source)
    dates = _daily_date_sequence(years)
    nsteps = length(dates)
    # One file per day — logging per-file here would be thousands of lines,
    # so just announce each field and how many days it covers up front.
    for name in fields
        @info "  loading $source $name ($nsteps daily files, $(first(dates)) to $(last(dates)))..."
    end
    rasters = _parallel_fetch(fields, dates) do name, d
        path = getraster(source, name; date = d, extras...)
        read(crop(Raster(path; lazy = true); to = area, touches = true))
    end
    layers = map(1:length(fields)) do fi
        per_day = @view rasters[fi, :]
        first_day = first(per_day)
        spatial_dims = dims(first_day)
        stacked = cat(map(parent, per_day)...; dims = 3)
        Raster(stacked, (spatial_dims..., Ti(1:nsteps)); crs = crs(first_day))
    end
    return NamedTuple{fields}(Tuple(layers))
end

# Full calendar-day sequence across `years` (Feb 29 included) — matches
# `_days_of_year`/`allocate_weather_buffers`'s leap-aware sizing.
function _daily_date_sequence(years)
    dates = Date[]
    for y in years
        for d in Date(y, 1, 1):Day(1):Date(y, 12, 31)
            push!(dates, d)
        end
    end
    return dates
end

function _load_layers(::HourlyZarrStore, source, fields::Tuple, area::Extent, years)
    # `getraster(source)` returns a `CachedCloudSource(url, cache_path)`.
    cloud_source = getraster(source)
    full_stack = RasterStack(cloud_source.url; source = Rasters.Zarrsource(), lazy = true)
    time_start = DateTime(first(years), 1, 1, 0)
    time_end = DateTime(last(years), 12, 31, 23)
    # `view` selects by dim name, so Zarr's storage order doesn't matter.
    layers = map(fields) do name
        long_name = layername(source, name)
        @info "  loading $source $name ($long_name)..."
        raw = getproperty(full_stack, Symbol(long_name))
        read(view(raw,
            X(area.X[1] .. area.X[2]),
            Y(area.Y[1] .. area.Y[2]),
            Ti(time_start .. time_end),
        ))
    end
    return NamedTuple{fields}(layers)
end

# Build one 3-D climatology raster (X, Y, Ti = 12 * nyears) from 12 monthly
# 2-D rasters via `path_for(month)`, tiled across `nyears`.
function _build_monthly_climatology(path_for, area::Extent, nyears::Int)
    monthly = map(1:12) do m
        path = path_for(m)
        read(crop(Raster(path; lazy = true); to = area, touches = true))
    end
    first_monthly = first(monthly)
    spatial_dims = dims(first_monthly)
    climatology = cat(map(parent, monthly)...; dims = 3)
    tiled = nyears == 1 ? climatology :
                                  repeat(climatology; outer = (1, 1, nyears))
    return Raster(tiled, (spatial_dims..., Ti(1:(12 * nyears))); crs = crs(first_monthly))
end

# Recursive unique-by-field reduction over a tuple of WeatherVariables.
# Type-stable: each step folds because `native_field(var)` is a constant
# Symbol pulled from the WeatherVariable's type parameter.
@inline _unique_native_fields(vars::Tuple) = _unique_native_fields((), vars)
@inline _unique_native_fields(acc::Tuple, ::Tuple{}) = acc
@inline function _unique_native_fields(acc::Tuple, vars::Tuple)
    f = native_field(first(vars))
    new_acc = _contains(acc, f) ? acc : (acc..., f)
    return _unique_native_fields(new_acc, Base.tail(vars))
end
@inline _contains(::Tuple{}, _) = false
@inline _contains(t::Tuple, f) = first(t) === f || _contains(Base.tail(t), f)

# ---------------------------------------------------------------------------
# Derivation registry
# ---------------------------------------------------------------------------

"""
    derive!(::Val{name}, buffers, ctx) -> nothing

Compute canonical variable `name` in place into `buffers.<name>`. Dispatched
by `Val(name)` so the driver can unroll over the derivation chain with one
compiled path per derivation.
"""
function derive! end

# Unified derivation chain, in topological order — every source's set of
# derivations is this list minus whatever it provides natively (e.g.
# TerraClimate skips wind_speed/vapour_pressure_deficit). Steps whose output
# isn't native but whose inputs are unpopulated just compute unused garbage;
# no runtime input-availability checks needed.
const _DEFAULT_DERIVATIONS = (
    Val(:reference_temperature_max), # ← lapse(maximum_temperature)
    Val(:reference_temperature_min), # ← lapse(minimum_temperature)
    Val(:mean_temperature), # ← (ref_max + ref_min) / 2
    Val(:wind_speed), # ← √(u_wind² + v_wind²)   — NCEP
    Val(:actual_vapour_pressure), # ← q·p / (0.622 + 0.378·q) — NCEP
    Val(:vapour_pressure_deficit), # ← e_s(mean_T) - actual_VP  — WorldClim, NCEP
    Val(:reference_humidity_max), # ← from VPD + ref_temp_min
    Val(:reference_humidity_min), # ← from VPD + ref_temp_max
    Val(:reference_wind_max), # ← wind_speed × shear_factor
    Val(:reference_wind_min), # ← ref_wind_max × 0.1
    Val(:cloud_cover), # ← from downward_shortwave_radiation
    Val(:cloud_min), # ← cloud_cover × 0.5 clamp
    Val(:cloud_max), # ← cloud_cover × 2.0 clamp
    Val(:deep_soil_temperature), # ← annual-mean of mean_temperature
)

# Hourly-mode chain — sources like ERA5 with hourly values directly. No
# min/max envelope; canonical hourly buffers *are* the env_hourly arrays.
const _DERIVATIONS_HOURLY = (
    Val(:wind_speed), # u/v → reference_wind_speed
    Val(:actual_vapour_pressure_from_dewpoint), # T_dew → actual_VP
    Val(:reference_humidity), # actual_VP + T → RH
    Val(:rainfall_daily_from_hourly), # hourly rainfall → daily total
    Val(:deep_soil_temperature_from_hourly), # annual mean of hourly T
)

# Hourly chain for sources with native RH (BARRA) instead of dewpoint
# (ERA5) — skips the dewpoint→VP→RH leg entirely.
const _DERIVATIONS_HOURLY_NATIVE_RH = (
    Val(:wind_speed),
    Val(:rainfall_daily_from_hourly),
    Val(:deep_soil_temperature_from_hourly),
)

# Sub-daily (6-hourly) chain — sources like NCEP whose values need
# interpolating/disaggregating to hourly.
const _DERIVATIONS_6H_TO_1H = (
    Val(:solar_geometry), # populate scratch.solar.out with hourly clear-sky
    Val(:wind_speed_6h), # u/v → scalar wind at 6h
    Val(:actual_vapour_pressure_6h), # specific_humidity + pressure → VP at 6h
    Val(:interpolate_met_6h_to_1h), # linear 6h→1h: T, wind, VP, pressure, LW
    Val(:apply_wind_shear),
    Val(:disaggregate_radiation_6h_to_1h), # solar-aware 6h→1h shortwave
    Val(:reference_humidity), # actual_VP + T → RH
    Val(:rainfall_daily_from_6h), # sum 4 × 6h blocks/day → daily total
    Val(:deep_soil_temperature_from_hourly), # annual mean
)

"""
    weather_derivations(::Type{<:RasterDataSource}) -> Tuple{Val, …}

Ordered derivation chain for this source. Defaults to
`_DEFAULT_DERIVATIONS`; hourly sources should override to
`_DERIVATIONS_HOURLY`.
"""
@inline weather_derivations(::Type) = _DEFAULT_DERIVATIONS

# ---------------------------------------------------------------------------
# Buffer allocation
# ---------------------------------------------------------------------------

"""
    allocate_weather_buffers(source, years) -> NamedTuple

Per-worker scratch for `assemble_weather!`: one unit-tagged buffer per
canonical variable, plus the env structs that share storage with them.

`years` is the real calendar-year span (not just a count) — leap-aware
buffer sizing needs it. Env-minmax struct type follows the source's
`temporal_resolution`.
"""
function allocate_weather_buffers(source::Type, years)
    resolution = temporal_resolution(source)
    days_of_year = _days_of_year(resolution, years)
    year_offsets = _year_offsets(resolution, years)
    _allocate_weather_buffers(resolution, source, days_of_year, year_offsets)
end

# Cumulative day/step boundaries per calendar year, e.g. `[0, 365, 730,
# 1096]` (leap years extend the relevant gap to 366). Built once per worker
# and read by reference thereafter, so consumers (`_annual_means!`,
# `derive!(:deep_soil_temperature_from_hourly)`) stay allocation-free.
function _year_offsets(resolution::TemporalResolution, years)
    lengths = _year_lengths(resolution, years)
    offsets = Vector{Int}(undef, length(lengths) + 1)
    offsets[1] = 0
    @inbounds for i in eachindex(lengths)
        offsets[i + 1] = offsets[i] + lengths[i]
    end
    return offsets
end

function _allocate_weather_buffers(resolution::Union{MonthlyResolution, DailyResolution},
                                   ::Any, days_of_year::AbstractVector{Int},
                                   year_offsets::Vector{Int})
    nsteps = length(days_of_year)
    zeros_float(n) = zeros(Float64, n)

    # Canonical buffers — each carries the unit appropriate for its quantity.
    # Some arrays are also held inside the env structs below; writes via the
    # canonical name show up in the struct (same array).
    maximum_temperature = zeros(typeof(0.0u"K"), nsteps)
    minimum_temperature = zeros(typeof(0.0u"K"), nsteps)
    mean_temperature = zeros(typeof(0.0u"K"), nsteps)
    reference_temperature_min = zeros(typeof(0.0u"K"), nsteps)
    reference_temperature_max = zeros(typeof(0.0u"K"), nsteps)
    wind_speed = zeros(typeof(0.0u"m/s"), nsteps)
    u_wind = zeros(typeof(0.0u"m/s"), nsteps)
    v_wind = zeros(typeof(0.0u"m/s"), nsteps)
    reference_wind_min = zeros(typeof(0.0u"m/s"), nsteps)
    reference_wind_max = zeros(typeof(0.0u"m/s"), nsteps)
    vapour_pressure_deficit = zeros(typeof(0.0u"kPa"), nsteps)
    actual_vapour_pressure = zeros(typeof(0.0u"kPa"), nsteps)
    specific_humidity = zeros_float(nsteps)
    surface_pressure = zeros(typeof(0.0u"Pa"), nsteps)
    reference_humidity_min = zeros_float(nsteps)
    reference_humidity_max = zeros_float(nsteps)
    downward_shortwave_radiation = zeros(typeof(0.0u"W/m^2"), nsteps)
    cloud_cover = zeros_float(nsteps)
    cloud_min = zeros_float(nsteps)
    cloud_max = zeros_float(nsteps)
    rainfall = zeros(typeof(0.0u"kg/m^2"), nsteps)
    deep_soil_temperature = zeros(typeof(0.0u"K"), nsteps)
    soil_moisture = zeros_float(nsteps)

    # Same field set for Monthly and Daily, so the kwargs are identical.
    environment_minmax = _minmax_env_type(resolution)(;
        reference_temperature_min, reference_temperature_max,
        reference_wind_min, reference_wind_max,
        reference_humidity_min, reference_humidity_max,
        cloud_min, cloud_max,
        minima_times = (temp = 0, wind = 0, humidity = 1, cloud = 1),
        maxima_times = (temp = 1, wind = 1, humidity = 0, cloud = 0),
    )
    # DailyTimeseries is the "per-step" struct regardless of resolution —
    # one entry per main timestep (month for monthly, day for daily).
    environment_daily = DailyTimeseries(;
        shade = zeros_float(nsteps),
        soil_wetness = zeros_float(nsteps),
        surface_emissivity = fill(0.95, nsteps),
        cloud_emissivity = fill(0.95, nsteps),
        rainfall, deep_soil_temperature,
        leaf_area_index = fill(0.1, nsteps),
    )
    # Pressure is constant per pixel; this is just a stub, rebuilt as a
    # `Fill` and swapped in per call. HourlyTimeseries is always 24× nsteps.
    environment_hourly = HourlyTimeseries(;
        pressure = Fill(atmospheric_pressure(0.0u"m"), nsteps * 24),
        reference_temperature = nothing, reference_humidity = nothing,
        reference_wind_speed = nothing, global_radiation = nothing,
        longwave_radiation = nothing, cloud_cover = nothing,
        rainfall = nothing, zenith_angle = nothing,
    )

    return (;
        maximum_temperature, minimum_temperature, mean_temperature,
        reference_temperature_min, reference_temperature_max,
        wind_speed, u_wind, v_wind, reference_wind_min, reference_wind_max,
        vapour_pressure_deficit, actual_vapour_pressure,
        specific_humidity, surface_pressure,
        reference_humidity_min, reference_humidity_max,
        downward_shortwave_radiation, cloud_cover, cloud_min, cloud_max,
        rainfall, deep_soil_temperature, soil_moisture,
        days_of_year, year_offsets,
        environment_minmax, environment_daily, environment_hourly,
    )
end

# Hourly mode: env_minmax = nothing, env_hourly arrays ARE the canonical
# buffers (shared storage), env_daily arrays sit at 1/24 the resolution
# and are filled by aggregating the hourly canonical buffers.
function _allocate_weather_buffers(::HourlyResolution, ::Any, days_of_year::AbstractVector{Int},
                                   year_offsets::Vector{Int})
    ndays  = length(days_of_year)    # one per calendar day (leap days included)
    nhours = ndays * 24
    zeros_float(n) = zeros(Float64, n)

    # Hourly canonical buffers — shared with `environment_hourly`.
    reference_temperature = zeros(typeof(0.0u"K"), nhours)
    reference_humidity = zeros_float(nhours)
    wind_speed = zeros(typeof(0.0u"m/s"), nhours)
    pressure = zeros(typeof(0.0u"Pa"), nhours)
    cloud_cover = zeros_float(nhours)
    global_radiation = zeros(typeof(0.0u"W/m^2"), nhours)
    longwave_radiation = zeros(typeof(0.0u"W/m^2"), nhours)
    rainfall = zeros(typeof(0.0u"kg/m^2"), nhours)
    zenith_angle = zeros(typeof(0.0u"°"), nhours)

    # Hourly intermediates feeding the small derivation chain.
    u_wind = zeros(typeof(0.0u"m/s"), nhours)
    v_wind = zeros(typeof(0.0u"m/s"), nhours)
    dewpoint_temperature = zeros(typeof(0.0u"K"), nhours)
    actual_vapour_pressure = zeros(typeof(0.0u"kPa"), nhours)
    # Mean sea level pressure — sources without native surface pressure
    # (BARRA's `:psl`) feed this into `derive!(:pressure_from_sea_level)`.
    sea_level_pressure = zeros(typeof(0.0u"Pa"), nhours)

    # `environment_daily` lives at one entry per unique day — aggregated from hourly.
    rainfall_daily = zeros(typeof(0.0u"kg/m^2"), ndays)
    deep_soil_temperature = zeros(typeof(0.0u"K"), ndays)
    soil_moisture = zeros_float(ndays)

    environment_daily = DailyTimeseries(;
        shade = zeros_float(ndays),
        soil_wetness = zeros_float(ndays),
        surface_emissivity = fill(0.95, ndays),
        cloud_emissivity = fill(0.95, ndays),
        rainfall = rainfall_daily,
        deep_soil_temperature,
        leaf_area_index = fill(0.1, ndays),
    )
    # All HourlyTimeseries fields populated — solver dispatches on the
    # `environment_minmax::Nothing` method and reads everything from here.
    environment_hourly = HourlyTimeseries(;
        pressure, reference_temperature, reference_humidity,
        reference_wind_speed = wind_speed,
        global_radiation, longwave_radiation, cloud_cover,
        rainfall, zenith_angle,
    )

    return (;
        reference_temperature, reference_humidity, wind_speed,
        u_wind, v_wind, pressure, cloud_cover,
        global_radiation, longwave_radiation,
        rainfall, rainfall_daily, zenith_angle,
        dewpoint_temperature, actual_vapour_pressure, sea_level_pressure,
        deep_soil_temperature, soil_moisture,
        days_of_year, year_offsets,
        environment_minmax = nothing,
        environment_daily, environment_hourly,
    )
end

# Sub-daily (6-hourly) mode: `*_6h` staging buffers (4/day, written by
# _read_native!) feed `_DERIVATIONS_6H_TO_1H` into hourly output buffers
# (24/day, shared with environment_hourly as in HourlyResolution), which in
# turn feed daily aggregates shared with environment_daily.
# TODO: staging field names differ from HourlyResolution's (u/v_wind_6h,
# specific_humidity_6h/surface_pressure_6h vs bare u/v_wind/dewpoint), so
# this can't just extend the hourly allocator — address in the
# SubDailyResolution{N} refactor.
function _allocate_weather_buffers(::SixHourlyResolution, _source,
                                   days_of_year::AbstractVector{Int},
                                   year_offsets::Vector{Int})
    ndays  = length(days_of_year)   # one per calendar day (leap days included)
    n6h    = ndays * 4              # native 6h staging steps
    nhours = ndays * 24             # hourly output
    zeros_float(n) = zeros(Float64, n)

    # 6h staging buffers (written by _read_native! via WeatherVariable declarations)
    reference_temperature_6h          = zeros(typeof(0.0u"K"),     n6h)
    u_wind_6h                         = zeros(typeof(0.0u"m/s"),   n6h)
    v_wind_6h                         = zeros(typeof(0.0u"m/s"),   n6h)
    specific_humidity_6h              = zeros_float(n6h)
    surface_pressure_6h               = zeros(typeof(0.0u"Pa"),    n6h)
    downward_shortwave_radiation_6h   = zeros(typeof(0.0u"W/m^2"), n6h)
    longwave_radiation_6h             = zeros(typeof(0.0u"W/m^2"), n6h)
    rainfall_6h                       = zeros(typeof(0.0u"kg/m^2"), n6h)
    # Derived at 6h resolution before interpolation
    wind_speed_6h                     = zeros(typeof(0.0u"m/s"),   n6h)
    actual_vapour_pressure_6h         = zeros(typeof(0.0u"kPa"),   n6h)

    # Hourly canonical output buffers — shared with environment_hourly
    reference_temperature = zeros(typeof(0.0u"K"),     nhours)
    reference_humidity    = zeros_float(nhours)
    wind_speed            = zeros(typeof(0.0u"m/s"),   nhours)
    pressure              = zeros(typeof(0.0u"Pa"),    nhours)
    cloud_cover           = zeros_float(nhours)
    global_radiation      = zeros(typeof(0.0u"W/m^2"), nhours)
    longwave_radiation    = zeros(typeof(0.0u"W/m^2"), nhours)
    rainfall              = zeros(typeof(0.0u"kg/m^2"), nhours)
    zenith_angle          = zeros(typeof(0.0u"°"),     nhours)
    actual_vapour_pressure = zeros(typeof(0.0u"kPa"),  nhours)

    # Daily aggregate buffers — shared with environment_daily
    rainfall_daily        = zeros(typeof(0.0u"kg/m^2"), ndays)
    deep_soil_temperature = zeros(typeof(0.0u"K"),       ndays)
    soil_moisture         = zeros_float(ndays)

    environment_daily = DailyTimeseries(;
        shade               = zeros_float(ndays),
        soil_wetness        = zeros_float(ndays),
        surface_emissivity  = fill(0.95, ndays),
        cloud_emissivity    = fill(0.95, ndays),
        rainfall            = rainfall_daily,
        deep_soil_temperature,
        leaf_area_index     = fill(0.1, ndays),
    )
    environment_hourly = HourlyTimeseries(;
        pressure, reference_temperature, reference_humidity,
        reference_wind_speed = wind_speed,
        global_radiation, longwave_radiation, cloud_cover,
        rainfall, zenith_angle,
    )

    return (;
        # 6h staging
        reference_temperature_6h, u_wind_6h, v_wind_6h,
        specific_humidity_6h, surface_pressure_6h,
        downward_shortwave_radiation_6h, longwave_radiation_6h, rainfall_6h,
        wind_speed_6h, actual_vapour_pressure_6h,
        # hourly output
        reference_temperature, reference_humidity, wind_speed,
        pressure, cloud_cover, global_radiation, longwave_radiation,
        rainfall, zenith_angle, actual_vapour_pressure,
        # daily
        rainfall_daily, deep_soil_temperature, soil_moisture,
        days_of_year, year_offsets,
        environment_minmax = nothing,
        environment_daily, environment_hourly,
    )
end

# ---------------------------------------------------------------------------
# Generic assembler
# ---------------------------------------------------------------------------

"""
    assemble_weather!(scratch, weather, source, site, I; kwargs...) -> NamedTuple

Build `(environment_minmax, environment_daily, environment_hourly)` for one
pixel of `weather` (a `RasterStack` loaded by `_load_weather(source, …)`).
Works in place into the scratch buffers held by `scratch.weather`.

  1. Read `source`'s native canonical variables (per `weather_variables`)
     into the canonical buffers, applying native-unit tagging and any
     source-specific transform.
  2. Run `weather_derivations(source)`, skipping any output that the
     source provides natively.
  3. Compute site-only quantities (pressure) and swap a fresh pressure
     `Fill` into the pre-built `environment_hourly` struct.
"""
function assemble_weather!(
    scratch, weather,
    source::Type,
    site, I::Tuple;
    grid_elevation = site.elevation,
    vapour_pressure_method = GoffGratch(),
    lapse_rate_model::LapseRate = EnvironmentalLapseRate(),
    canonical_overrides::NamedTuple = (;),
    wind_reference_height = 2.0u"m",
)
    buffers = scratch.weather
    atm_pressure = atmospheric_pressure(site.elevation)
    variables = weather_variables(source)

    ctx = (;
        site, grid_elevation, lapse_rate_model, vapour_pressure_method,
        atmospheric_pressure = atm_pressure, scratch,
        wind_reference_height,
    )

    _read_native!(buffers, weather, variables, I)
    _run_derivations!(buffers, ctx, source, variables, canonical_overrides)
    _apply_canonical_overrides!(buffers, canonical_overrides, I)

    environment_hourly = _finalize_environment_hourly(buffers, variables, atm_pressure)

    return (; buffers.environment_minmax, buffers.environment_daily, environment_hourly)
end

# Monthly/daily mode: pressure is constant per pixel, so rebuild the
# env_hourly `Fill` from site elevation and swap it in. Hourly/sub-daily
# mode already has a populated pressure array (native or interpolated from
# 6h staging) that must not be overwritten — `environment_minmax ===
# nothing` unambiguously identifies that case (set only by those allocators).
@inline function _finalize_environment_hourly(buffers, variables, atm_pressure)
    buffers.environment_minmax === nothing && return buffers.environment_hourly
    _is_native(Val(:pressure), variables)  && return buffers.environment_hourly
    pressure = Fill(atm_pressure, length(buffers.environment_hourly.pressure))
    return setproperties(buffers.environment_hourly, (; pressure))
end

# ---------------------------------------------------------------------------
# Native-variable read step
# ---------------------------------------------------------------------------

@inline function _read_native!(buffers, weather, variables::Tuple, I)
    unrolled_map(variables) do var
        _read_one_variable!(buffers, weather, var, I)
        nothing
    end
    return nothing
end

@inline function _read_one_variable!(buffers, weather,
                                     var::WeatherVariable{Name, Field}, I) where {Name, Field}
    target = getproperty(buffers, Name)
    layer = getproperty(weather, Field)
    transform = var.transform
    unit = var.unit
    # Splat the spatial dim tuple from `DimIndices`; `Ti(k)` selects the
    # k-th timestep regardless of how `layer` stores its dims.
    @inbounds for k in eachindex(target)
        target[k] = transform(layer[I..., Ti(k)]) * unit
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Derivation driver
# ---------------------------------------------------------------------------

@inline function _run_derivations!(buffers, ctx, source, variables::Tuple,
                                   overrides::NamedTuple = (;))
    unrolled_map(weather_derivations(source)) do v
        skip = _is_native(v, variables) || _is_override(v, overrides)
        skip || derive!(v, buffers, ctx)
        nothing
    end
    return nothing
end

# True iff the canonical variable named `N` has a user-supplied override
# raster — checked at compile time via the NamedTuple's `K` parameter so
# the result folds into the derivation skip.
@inline _is_override(::Val{N}, ::NamedTuple{K}) where {N, K} = N in K

# Apply per-canonical-variable user overrides. Runs after derivations so
# overrides win even if the source provides the variable natively.
@inline function _apply_canonical_overrides!(buffers, overrides::NamedTuple{K}, I) where K
    unrolled_map(K) do name
        _apply_one_override!(buffers, name, getproperty(overrides, name), I)
        nothing
    end
    return nothing
end
@inline _apply_canonical_overrides!(_, ::NamedTuple{()}, _) = nothing

@inline function _apply_one_override!(buffers, name::Symbol, raster, I)
    target = getproperty(buffers, name)
    _copy_override!(target, raster, I)
    return nothing
end

# Dispatches on Ti presence: constant-in-time override is a spatial-only
# Raster (rank == length(I)); time-varying has a trailing Ti dim (rank ==
# length(I) + 1). `hasdim` is a compile-time check, so this folds away.
@inline function _copy_override!(target::AbstractVector, raster, I)
    @inbounds if hasdim(raster, Ti)
        for k in eachindex(target)
            target[k] = raster[I..., Ti(k)]
        end
    else
        value = raster[I...]
        for k in eachindex(target)
            target[k] = value
        end
    end
    return nothing
end

# Recursive check over the variables tuple — each step dispatches on the
# WeatherVariable's Name type param so the comparison folds at compile time.
@inline _is_native(::Val, ::Tuple{}) = false
@inline _is_native(v::Val{N}, vars::Tuple) where {N} =
    _name_matches(v, first(vars)) || _is_native(v, Base.tail(vars))

@inline _name_matches(::Val{N}, ::WeatherVariable{N}) where {N} = true
@inline _name_matches(::Val, ::WeatherVariable) = false

# ---------------------------------------------------------------------------
# Derivations
# ---------------------------------------------------------------------------

function derive!(::Val{:reference_temperature_max}, buffers, ctx)
    _lapse_correct!(buffers.reference_temperature_max,
                    buffers.maximum_temperature, ctx)
end

function derive!(::Val{:reference_temperature_min}, buffers, ctx)
    _lapse_correct!(buffers.reference_temperature_min,
                    buffers.minimum_temperature, ctx)
end

@inline function _lapse_correct!(out, src, ctx)
    (; site, grid_elevation, lapse_rate_model) = ctx
    Δz = site.elevation - grid_elevation
    if iszero(Δz)
        @. out = src
    else
        lr = lapse_rate(lapse_rate_model)   # scalar K/m — extract before broadcast
        @. out = src - lr * Δz
    end
    return nothing
end

function derive!(::Val{:mean_temperature}, buffers, ctx)
    @. buffers.mean_temperature =
        (buffers.reference_temperature_max + buffers.reference_temperature_min) / 2
    return nothing
end

# q · p / (0.622 + 0.378 · q), for sources (NCEP) providing specific
# humidity + surface pressure rather than vapour pressure directly.
function derive!(::Val{:actual_vapour_pressure}, buffers, ctx)
    @inbounds for k in eachindex(buffers.actual_vapour_pressure)
        q = buffers.specific_humidity[k]
        p = buffers.surface_pressure[k]
        buffers.actual_vapour_pressure[k] = q * p / (0.622 + 0.378 * q)
    end
    return nothing
end

# VPD = e_s(mean_T) - actual_VP, for sources with actual vapour pressure
# (WorldClim's `:vapr`, or NCEP via the SH→VP derivation above).
function derive!(::Val{:vapour_pressure_deficit}, buffers, ctx)
    method = ctx.vapour_pressure_method
    @inbounds for k in eachindex(buffers.vapour_pressure_deficit)
        saturation = vapour_pressure(method, buffers.mean_temperature[k])
        buffers.vapour_pressure_deficit[k] = saturation - buffers.actual_vapour_pressure[k]
    end
    return nothing
end

function derive!(::Val{:reference_humidity_max}, buffers, ctx)
    # NicheMapR micro_terra.R pairing: RH_max with T_min.
    _relative_humidity_from_vpd!(buffers.reference_humidity_max,
        buffers.vapour_pressure_deficit, buffers.mean_temperature,
        buffers.reference_temperature_min, ctx.vapour_pressure_method)
    return nothing
end

function derive!(::Val{:reference_humidity_min}, buffers, ctx)
    # NicheMapR micro_terra.R pairing: RH_min with T_max.
    _relative_humidity_from_vpd!(buffers.reference_humidity_min,
        buffers.vapour_pressure_deficit, buffers.mean_temperature,
        buffers.reference_temperature_max, ctx.vapour_pressure_method)
    return nothing
end

# Vector magnitude of U/V wind components — used by sources (NCEP, ERA5)
# that provide horizontal wind as east-west and north-south components
# rather than as a scalar speed.
function derive!(::Val{:wind_speed}, buffers, ctx)
    @. buffers.wind_speed = sqrt(buffers.u_wind^2 + buffers.v_wind^2)
    return nothing
end

# Power-law height correction: scale wind from `z_src` to `z_ref` with
# neutral-stability shear exponent α = 0.15.
_wind_height_correction(z_ref, z_src = 10.0u"m", α = 0.15) = (z_ref / z_src)^α

function derive!(::Val{:reference_wind_max}, buffers, ctx)
    shear = _wind_height_correction(ctx.wind_reference_height)
    @. buffers.reference_wind_max = buffers.wind_speed * shear
    return nothing
end

function derive!(::Val{:reference_wind_min}, buffers, ctx)
    @. buffers.reference_wind_min = buffers.reference_wind_max * 0.1
    return nothing
end

function derive!(::Val{:cloud_cover}, buffers, ctx)
    (; scratch, site) = ctx
    # `flat_terrain_template` is constant except for elevation/pressure/lat/lon;
    # swap those four in via `setproperties` (zero-cost).
    flat_terrain = setproperties(scratch.cloud_constants.flat_terrain_template,
        (; site.elevation,
           atmospheric_pressure = ctx.atmospheric_pressure,
           site.latitude, site.longitude),
    )
    cloud_from_solar_radiation!(buffers.cloud_cover, scratch.solar.out,
        buffers.downward_shortwave_radiation, flat_terrain, buffers.days_of_year,
        scratch.solar.buffers;
        solar_model = scratch.cloud_constants.solar_model,
        hours = scratch.cloud_constants.hours,
    )
    return nothing
end

function derive!(::Val{:cloud_min}, buffers, ctx)
    @. buffers.cloud_min = clamp(buffers.cloud_cover * 0.5, 0.0, 1.0)
    return nothing
end

function derive!(::Val{:cloud_max}, buffers, ctx)
    @. buffers.cloud_max = clamp(buffers.cloud_cover * 2.0, 0.0, 1.0)
    return nothing
end

function derive!(::Val{:deep_soil_temperature}, buffers, ctx)
    _annual_means!(buffers.deep_soil_temperature,
                   buffers.mean_temperature, buffers.year_offsets)
    return nothing
end

# ---------------------------------------------------------------------------
# Hourly-mode derivations
# ---------------------------------------------------------------------------

# At the dewpoint the air is saturated, so actual VP = saturation VP at
# the dewpoint (ERA5's 2 m dewpoint temperature).
function derive!(::Val{:actual_vapour_pressure_from_dewpoint}, buffers, ctx)
    method = ctx.vapour_pressure_method
    @inbounds for k in eachindex(buffers.actual_vapour_pressure)
        buffers.actual_vapour_pressure[k] =
            vapour_pressure(method, buffers.dewpoint_temperature[k])
    end
    return nothing
end

# Barometric formula at site elevation, using each hour's mean sea level
# pressure (BARRA's `:psl`) as reference rather than the standard-atmosphere
# default — BARRA's grid elevation differs from the DEM-derived site one.
function derive!(::Val{:pressure_from_sea_level}, buffers, ctx)
    elevation = ctx.site.elevation
    @inbounds for k in eachindex(buffers.pressure)
        buffers.pressure[k] =
            atmospheric_pressure(elevation; reference_pressure = buffers.sea_level_pressure[k])
    end
    return nothing
end

# Single-value hourly relative humidity:
#   RH = actual_VP / saturation_VP(reference_temperature)
# Output clamped to [0, 1].
function derive!(::Val{:reference_humidity}, buffers, ctx)
    method = ctx.vapour_pressure_method
    @inbounds for k in eachindex(buffers.reference_humidity)
        saturation = vapour_pressure(method, buffers.reference_temperature[k])
        if saturation <= zero(saturation)
            buffers.reference_humidity[k] = 1.0
        else
            buffers.reference_humidity[k] =
                clamp(buffers.actual_vapour_pressure[k] / saturation, 0.0, 1.0)
        end
    end
    return nothing
end

# Aggregate 24 hourly rainfall values into one daily total for env_daily.
# TODO: propagate hourly rainfall into Microclimate.jl instead — the solver
# should consume `env_hourly.rainfall` directly; pre-aggregating here loses
# sub-daily intensity info that affects infiltration/runoff.
function derive!(::Val{:rainfall_daily_from_hourly}, buffers, ctx)
    ndays = length(buffers.rainfall_daily)
    @inbounds for d in 1:ndays
        h0 = (d - 1) * 24
        s = zero(eltype(buffers.rainfall))
        for h in 1:24
            s += buffers.rainfall[h0 + h]
        end
        buffers.rainfall_daily[d] = s
    end
    return nothing
end

# Annual mean of hourly reference_temperature, broadcast across each year's
# `deep_soil_temperature` entries. `year_offsets` gives variable-length
# (leap-aware) year blocks without allocating inside this per-pixel derivation.
function derive!(::Val{:deep_soil_temperature_from_hourly}, buffers, ctx)
    year_offsets = buffers.year_offsets
    ndays  = length(buffers.deep_soil_temperature)
    nhours = length(buffers.reference_temperature)
    nyears = length(year_offsets) - 1
    if nyears <= 0 || year_offsets[end] != ndays
        # Sub-annual or mismatched run: fill with mean of all available hours.
        run_mean = sum(buffers.reference_temperature) / nhours
        fill!(buffers.deep_soil_temperature, run_mean)
        return nothing
    end
    @inbounds for y in 1:nyears
        d_lo, d_hi = year_offsets[y] + 1, year_offsets[y + 1]
        h_lo, h_hi = (d_lo - 1) * 24 + 1, d_hi * 24
        annual_sum = zero(eltype(buffers.reference_temperature))
        for h in h_lo:h_hi
            annual_sum += buffers.reference_temperature[h]
        end
        annual_mean = annual_sum / (h_hi - h_lo + 1)
        for d in d_lo:d_hi
            buffers.deep_soil_temperature[d] = annual_mean
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Sub-daily (6-hourly) derivations
# ---------------------------------------------------------------------------

# Populates scratch.solar.out with hourly clear-sky global_horizontal for
# all days. Must run first in _DERIVATIONS_6H_TO_1H.
function derive!(::Val{:solar_geometry}, buffers, ctx)
    (; scratch, site) = ctx
    flat_terrain = setproperties(scratch.cloud_constants.flat_terrain_template,
        (; site.elevation,
           atmospheric_pressure = ctx.atmospheric_pressure,
           site.latitude, site.longitude))
    solar_radiation!(scratch.solar.out, scratch.solar.buffers,
        scratch.cloud_constants.solar_model;
        solar_terrain = flat_terrain,
        days  = buffers.days_of_year,
        hours = scratch.cloud_constants.hours)
    return nothing
end

# Scalar wind speed from U/V components at 6h resolution (at source height, 10 m).
function derive!(::Val{:wind_speed_6h}, buffers, _ctx)
    @. buffers.wind_speed_6h = sqrt(buffers.u_wind_6h^2 + buffers.v_wind_6h^2)
    return nothing
end

# Applied post-interpolation so wind_speed_6h stays a single concern.
function derive!(::Val{:apply_wind_shear}, buffers, ctx)
    shear = _wind_height_correction(ctx.wind_reference_height)
    @. buffers.wind_speed *= shear
    return nothing
end

# Actual vapour pressure from specific humidity + surface pressure at 6h resolution.
function derive!(::Val{:actual_vapour_pressure_6h}, buffers, _ctx)
    @inbounds for k in eachindex(buffers.actual_vapour_pressure_6h)
        q = buffers.specific_humidity_6h[k]
        p = buffers.surface_pressure_6h[k]
        buffers.actual_vapour_pressure_6h[k] = q * p / (0.622 + 0.378 * q)
    end
    return nothing
end

# 4 6-hourly values/day → 24 hourly, linearly interpolated between block
# midpoints (hours 3, 9, 15, 21); edges clamp to the boundary value. `ndays`
# comes from the buffer length, so this is leap-year-safe automatically.
function _interp_6h_to_1h!(out::AbstractVector, src::AbstractVector)
    ndays = length(src) ÷ 4
    @inbounds for d in 1:ndays
        k0 = (d - 1) * 4   # 0-based 6h offset
        v1 = src[k0 + 1];  v2 = src[k0 + 2]
        v3 = src[k0 + 3];  v4 = src[k0 + 4]
        # Previous/next boundary values for interp across day boundaries
        vp = d > 1    ? src[(d-2)*4 + 4] : v1
        vn = d < ndays ? src[d*4 + 1]    : v4
        h0 = (d - 1) * 24   # 0-based hour offset
        for h in 0:23
            hf = h + 0.5   # midpoint of each clock hour
            local op::typeof(v1)
            if hf < 3.0
                t = (hf + 3.0) / 6.0   # 0 at h=-3 (prev centroid 21), 1 at h=3
                op = vp * (1.0 - t) + v1 * t
            elseif hf < 9.0
                t = (hf - 3.0) / 6.0
                op = v1 * (1.0 - t) + v2 * t
            elseif hf < 15.0
                t = (hf - 9.0) / 6.0
                op = v2 * (1.0 - t) + v3 * t
            elseif hf < 21.0
                t = (hf - 15.0) / 6.0
                op = v3 * (1.0 - t) + v4 * t
            else
                t = (hf - 21.0) / 6.0
                op = v4 * (1.0 - t) + vn * t
            end
            out[h0 + h + 1] = op
        end
    end
    return out
end

# T, wind, humidity, pressure, longwave: 6h to hourly via _interp_6h_to_1h!.
function derive!(::Val{:interpolate_met_6h_to_1h}, buffers, _ctx)
    _interp_6h_to_1h!(buffers.reference_temperature,    buffers.reference_temperature_6h)
    _interp_6h_to_1h!(buffers.wind_speed,               buffers.wind_speed_6h)
    _interp_6h_to_1h!(buffers.actual_vapour_pressure,   buffers.actual_vapour_pressure_6h)
    _interp_6h_to_1h!(buffers.pressure,                 buffers.surface_pressure_6h)
    _interp_6h_to_1h!(buffers.longwave_radiation,       buffers.longwave_radiation_6h)
    return nothing
end

# Solar-aware shortwave disaggregation, 6h to hourly (NicheMapR/microclima
# approach): per 6h block, opacity = observed / clear-sky mean (from
# derive!(:solar_geometry)), applied constant within the block against the
# hourly clear-sky. Zero clear-sky mean (polar night) → zero radiation.
function derive!(::Val{:disaggregate_radiation_6h_to_1h}, buffers, ctx)
    gh = ctx.scratch.solar.out.global_horizontal   # 24×ndays hourly clear-sky
    ndays = length(buffers.days_of_year)           # 365 × nyears
    @inbounds for d in 1:ndays
        h0 = (d - 1) * 24   # 0-based index into gh (24 per day)
        for b in 0:3         # 4 six-hour blocks per day
            k6h = (d - 1) * 4 + b + 1   # 1-based 6h index
            # Mean clear-sky over this 6h block (hours b*6+1 … b*6+6)
            cs_sum = zero(eltype(gh))
            for j in 1:6
                cs_sum += gh[h0 + b*6 + j]
            end
            cs_mean = cs_sum / 6
            obs = buffers.downward_shortwave_radiation_6h[k6h]
            opacity = cs_mean <= zero(cs_mean) ? 0.0 :
                clamp(ustrip(u"W/m^2", obs) / ustrip(u"W/m^2", cs_mean), 0.0, 1.0)
            # Apply constant opacity to each hour in the block
            for j in 1:6
                hi = h0 + b*6 + j   # 1-based index into hourly arrays
                buffers.global_radiation[hi] = opacity * gh[h0 + b*6 + j]
            end
        end
    end
    return nothing
end

# Sum 4 six-hourly rainfall accumulations per day into daily totals.
function derive!(::Val{:rainfall_daily_from_6h}, buffers, _ctx)
    ndays = length(buffers.rainfall_daily)
    @inbounds for d in 1:ndays
        k0 = (d - 1) * 4
        buffers.rainfall_daily[d] = buffers.rainfall_6h[k0 + 1] +
                                    buffers.rainfall_6h[k0 + 2] +
                                    buffers.rainfall_6h[k0 + 3] +
                                    buffers.rainfall_6h[k0 + 4]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Helpers reused by derivations
# ---------------------------------------------------------------------------

"""
    _relative_humidity_from_vpd!(out, vapour_pressure_deficit,
                                 mean_temperature, reference_temperature, method)

NicheMapR micro_terra.R: actual_VP = saturation_VP(mean_T) − VPD;
RH = actual_VP / saturation_VP(reference_T), clamped to [0, 1]. Pass
`reference_temperature = minimum_temperature` for RH_max, `= maximum_temperature`
for RH_min.
"""
function _relative_humidity_from_vpd!(
    out::AbstractVector,
    vapour_pressure_deficit::AbstractVector,
    mean_temperature::AbstractVector,
    reference_temperature::AbstractVector,
    method,
)
    @inbounds for k in eachindex(out, vapour_pressure_deficit,
                                 mean_temperature, reference_temperature)
        saturation_at_mean = vapour_pressure(method, mean_temperature[k])
        actual = saturation_at_mean - vapour_pressure_deficit[k]
        if actual <= zero(actual)
            out[k] = 0.0
            continue
        end
        saturation_at_reference = vapour_pressure(method, reference_temperature[k])
        out[k] = saturation_at_reference <= zero(saturation_at_reference) ? 1.0 :
                 clamp(actual / saturation_at_reference, 0.0, 1.0)
    end
    return out
end

"""
    _annual_means!(out, values, year_offsets)

Replaces every entry in each calendar year's block (per `year_offsets`,
see `_year_offsets`) with that year's mean. Builds `deep_soil_temperature`
for monthly/daily/hourly-aggregated resolutions; leap-aware and
allocation-free since `year_offsets` is precomputed once per worker.
"""
function _annual_means!(out::AbstractVector,
                        values::AbstractVector,
                        year_offsets::Vector{Int})
    nyears = length(year_offsets) - 1
    if nyears <= 0 || year_offsets[end] != length(values)
        # Sub-annual or mismatched run: fill with mean of all available steps.
        run_mean = sum(values) / length(values)
        fill!(out, run_mean)
        return out
    end
    @inbounds for y in 1:nyears
        lo, hi = year_offsets[y] + 1, year_offsets[y + 1]
        annual_sum = zero(eltype(values))
        for k in lo:hi
            annual_sum += values[k]
        end
        annual_mean = annual_sum / (hi - lo + 1)
        for k in lo:hi
            out[k] = annual_mean
        end
    end
    return out
end
