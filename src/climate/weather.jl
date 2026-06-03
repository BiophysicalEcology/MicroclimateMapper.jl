# Generic per-pixel weather assembler.
#
# A weather source declares
#   * `weather_loader(source)`    — how to read the raster files (yearly
#                                   time-series, 12-file monthly climatology,
#                                   12-file future climatology, one multi-band
#                                   future file, …).
#   * `weather_variables(source)` — which canonical physical variables each
#                                   of the source's raster layers provides
#                                   (with native unit and optional transform).
#
# The generic `_load_weather` and `assemble_weather!` then do the rest: load
# the rasters via the declared loader, optionally merge fallback layers from
# a baseline source, read each native canonical variable into a unit-tagged
# buffer, and run a chain of derivation steps for every canonical variable
# the source doesn't natively provide.
#
# Adding a new source is just three small methods (loader, variables, and
# optionally fallback_source/fallback_layers/primary_layers for the cases
# where a source needs to fall back to a baseline for missing fields).

# ---------------------------------------------------------------------------
# Variable map: canonical → native (per source)
# ---------------------------------------------------------------------------

"""
    WeatherVariable{Name, Field, U, T}(unit, transform)

Declares that the source's raster layer `Field` provides the canonical
variable `Name`, with native `unit` (Unitful — auto-converted on
assignment) and optional `transform` for source-specific scaling that
isn't a unit conversion. Both `Name` and `Field` are encoded as type
parameters so a tuple of `WeatherVariable`s gives the compiler one
specialised path per entry.
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

Singleton trait describing the time step a source provides — currently
monthly (12 timesteps per year) or daily (365 timesteps per year). Drives
both the size of every per-pixel buffer and which env-minmax struct gets
built (`MonthlyMinMaxEnvironment` vs `DailyMinMaxEnvironment`).
"""
abstract type TemporalResolution end

"""
    MonthlyResolution

12 timesteps per simulated year. Used by TerraClimate, CHELSA, WorldClim.
Picks `MonthlyMinMaxEnvironment` and `Microclimate.DEFAULT_DAYS` (12 mid-month
days of year) for solar geometry.
"""
struct MonthlyResolution <: TemporalResolution end

"""
    DailyResolution

365 timesteps per simulated year (leap days dropped for now). Used by
daily sources such as NCEP, AWAP, ERA5 daily aggregates. Picks
`DailyMinMaxEnvironment` and `1:365` for solar geometry — the solver then
runs in consecutive-day mode (each day inherits soil state from the
previous, see `DailyMinMaxEnvironment` docstring in Microclimate).
"""
struct DailyResolution <: TemporalResolution end

"""
    HourlyResolution

8760 timesteps per simulated year (24 × 365; leap-day Feb 29 dropped to
keep year boundaries clean). Used by hourly sources such as ERA5. Skips
the min/max envelope entirely (`environment_minmax = nothing`) — hourly
values flow straight into `environment_hourly`. `environment_daily`
arrays remain at 365 entries per year and are filled by aggregating the
hourly canonical buffers.
"""
struct HourlyResolution <: TemporalResolution end

"""
    SixHourlyResolution

Native 6-hourly files (4 × 365 = 1460 steps/year) interpolated to 8760 hourly
output steps/year before the solver sees the data. Used by `NCEPHourly`. Like
`HourlyResolution`, skips the min/max envelope and runs in consecutive-day mode.

Future sub-daily sources (3-hourly MERRA-2, 12-hourly GFS, etc.) follow the
same pattern — declare their own `ThreeHourlyResolution` etc. with the matching
buffer/derivation-chain dispatches.
"""
struct SixHourlyResolution <: TemporalResolution end

"""
    temporal_resolution(::Type{<:RasterDataSource}) -> TemporalResolution

The time step a source's `_load_weather` produces. Default
`MonthlyResolution()`; daily/hourly sources must override.
"""
@inline temporal_resolution(::Type) = MonthlyResolution()

@inline steps_per_year(::MonthlyResolution) = 12
@inline steps_per_year(::DailyResolution) = 365
@inline steps_per_year(::HourlyResolution) = 8760
@inline steps_per_year(::SixHourlyResolution) = 8760  # output to solver is hourly

@inline _minmax_env_type(::MonthlyResolution) = MonthlyMinMaxEnvironment
@inline _minmax_env_type(::DailyResolution) = DailyMinMaxEnvironment

@inline _days_of_year(::MonthlyResolution, nyears) = repeat(Microclimate.DEFAULT_DAYS, nyears)
@inline _days_of_year(::DailyResolution, nyears) = repeat(1:365, nyears)
@inline _days_of_year(::HourlyResolution, nyears) = repeat(1:365, nyears)
@inline _days_of_year(::SixHourlyResolution, nyears) = repeat(1:365, nyears)

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

# Day-of-year on the 365-day calendar (Feb 29 is dropped, so March 1 = doy 60
# whether or not the year is a leap year).
function _doy_noleap(d::Date)
    doy = dayofyear(d)
    isleapyear(year(d)) && d >= Date(year(d), 3, 1) && (doy -= 1)
    return doy
end

_is_leapday(d::Date) = isleapyear(year(d)) && month(d) == 2 && day(d) == 29

# Normalise user-supplied dates to a sorted Vector{Date}, dropping Feb 29
# (consistent with the 365-day-per-year convention throughout).
_normalise_dates(d::Date) = _is_leapday(d) ? Date[] : [d]
function _normalise_dates(r::AbstractRange{Date})
    v = filter(!_is_leapday, collect(r))
    isempty(v) && error(
        "No valid simulation dates remain after dropping Feb 29 " *
        "(leap days are not yet supported).")
    return v
end

# Compute the 1-based positional Ti range into the full-years weather stack,
# and the day-of-year vector for the solver (one entry per unique day, or per
# month for monthly sources). `dates_vec` must already be normalised (sorted,
# no Feb 29).
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
    ti_start = (yi_start - 1) * 365 + _doy_noleap(start_d)
    ti_end   = (yi_end   - 1) * 365 + _doy_noleap(end_d)
    return ti_start, ti_end, _doy_noleap.(dates_vec)
end

function _ti_range_for_dates(::HourlyResolution, years, dates_vec::Vector{Date})
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = (yi_start - 1) * 8760 + (_doy_noleap(start_d) - 1) * 24 + 1
    ti_end   = (yi_end   - 1) * 8760 + _doy_noleap(end_d) * 24
    return ti_start, ti_end, _doy_noleap.(dates_vec)
end

# The loaded NCEPHourly stack is at native 6-hourly resolution (4 × 365 = 1460
# steps/year), so the Ti slice is in 6h units. The returned day-of-year vector
# is one entry per calendar day (matching the hourly solver output), as for
# `HourlyResolution`.
function _ti_range_for_dates(::SixHourlyResolution, years, dates_vec::Vector{Date})
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = (yi_start - 1) * 1460 + (_doy_noleap(start_d) - 1) * 4 + 1
    ti_end   = (yi_end   - 1) * 1460 + _doy_noleap(end_d) * 4
    return ti_start, ti_end, _doy_noleap.(dates_vec)
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

Singleton trait describing how a weather source's raster files are
organised. `_load_layers(::WeatherLoader, source, fields, area, years)`
turns a tuple of source-field symbols into a `NamedTuple` of 3-D rasters
`(X, Y, Ti = 12 * nyears)`.
"""
abstract type WeatherLoader end

"""
    YearlyTimeSeries

One file per year, each file containing all 12 months along a `Ti` axis.
Used by TerraClimate. `getraster(source, name; date = Date(year))`.
"""
struct YearlyTimeSeries <: WeatherLoader end

"""
    MonthlyClimatology

12 files per layer (one per month) representing a fixed climatology
(no time dimension in the data — tiled across the requested `years`).
Used by CHELSA{Climate} and WorldClim{Climate}.
`getraster(source, name; month)`.
"""
struct MonthlyClimatology <: WeatherLoader end

"""
    FutureMonthlyClimatology

12 files per layer for a specific future-period date. Like
`MonthlyClimatology` but with an extra `date` kwarg selecting the
projected window. Used by CHELSA{Future{Climate}}.
`getraster(source, name; date = future_date, month)`.
"""
struct FutureMonthlyClimatology <: WeatherLoader end

"""
    MultiBandFutureClimatology

One multi-band file per layer for a specific future-period date — the 12
months live inside a single GeoTIFF as separate bands. Used by
WorldClim{Future{Climate}}. `getraster(source, name; date = future_date)`.
"""
struct MultiBandFutureClimatology <: WeatherLoader end

"""
    DailyFiles

One file per day per layer — the source ships one 2-D raster per calendar
day. Used by daily station/grid products such as AWAP.
`getraster(source, name; date = day)`. Feb 29 in leap years is skipped
to match the fixed `nyears * 365` buffer size.
"""
struct DailyFiles <: WeatherLoader end

"""
    HourlyZarrStore

Single remote Zarr store containing every layer as a 3-D `(time, lat, lon)`
variable. Used by ERA5 via the ARCO-ERA5 cloud-optimised store. The
source's `getraster(source)` (no args) returns a `CachedCloudSource` whose
`url` is opened once via Rasters' ZarrDatasets extension and each layer
is read out by its long name (`RasterDataSources.layername(source, sym)`).
"""
struct HourlyZarrStore <: WeatherLoader end

"""
    weather_loader(::Type{<:RasterDataSource}) -> WeatherLoader

Singleton declaring how this source's files are organised. One of
`YearlyTimeSeries`, `MonthlyClimatology`, `FutureMonthlyClimatology`,
`MultiBandFutureClimatology`.
"""
function weather_loader end

# ---------------------------------------------------------------------------
# Per-source extension points (with defaults)
# ---------------------------------------------------------------------------

"""
    primary_layers(source) -> Tuple{Symbol, …}

Source-native raster layer names to load via `weather_loader(source)`.
Defaults to the unique set of `native_field`s drawn from
`weather_variables(source)`. Future-style sources whose primary file set
is a strict subset of their declared variables (the rest coming from a
baseline) should override this.
"""
@inline primary_layers(source) = _unique_native_fields(weather_variables(source))

"""
    fallback_source(source) -> source or Nothing

If `source` can't provide every canonical variable itself, declare the
baseline source to fall back on (e.g. CHELSA{Climate} for CHELSA Future).
Default `nothing` means no fallback.
"""
@inline fallback_source(::Type) = nothing

"""
    fallback_layers(source) -> Tuple{Symbol, …}

Source-field names to load from the `fallback_source` and merge into the
returned stack. Default empty.
"""
@inline fallback_layers(::Type) = ()

"""
    _extra_getraster_kwargs(source) -> NamedTuple

Extra keyword arguments to splat into every `getraster(source, name; …)`
call (in addition to the loader's own kwargs like `date`/`month`). NCEP
uses this to pass `dataset = "reanalysis"`. Default empty.
"""
@inline _extra_getraster_kwargs(::Type) = (;)

# ---------------------------------------------------------------------------
# Generic _load_weather
# ---------------------------------------------------------------------------

"""
    _post_load_stack!(source, stack, nyears) -> RasterStack

Source-specific post-processing hook called immediately after every layer of
`stack` has been loaded. Default is a no-op (returns `stack` unchanged).

Override for sources whose files contain sub-daily Ti that the declared
`temporal_resolution` does not expect — e.g. `NCEP{SurfaceGauss}` in daily
mode stores `dswrf` at 3-hourly resolution and needs it averaged to daily.
"""
_post_load_stack!(::Type, stack, _nyears) = stack

# Average a sub-daily Ti layer down to daily means.
# `factor` is the number of sub-daily steps per day (e.g. 8 for 3-hourly).
function _aggregate_ti_to_daily(layer::AbstractRaster, factor::Int)
    ndays = size(layer, Ti) ÷ factor
    # Reshape Ti into (factor, ndays) blocks and average along the first axis.
    xy_dims = dims(layer)[1:2]
    data = parent(layer)
    # Materialise to a plain Array so we can reshape.
    arr = Array(data)
    # Layout: (nx, ny, nsteps). Reshape last dim to (factor, ndays).
    nx, ny = size(arr, 1), size(arr, 2)
    reshaped = reshape(arr, nx, ny, factor, ndays)
    daily = dropdims(sum(reshaped; dims = 3); dims = 3) ./ factor
    return Raster(daily, (xy_dims..., Ti(1:ndays)); crs = crs(layer))
end

"""
    _load_weather(source, area, years) -> RasterStack

Loads every layer the source contributes to the canonical variable map.
For sources with a `fallback_source`, the primary layers come from the
source itself and the fallback layers come from the baseline — both are
merged into one `RasterStack`. Dispatch is via `weather_loader(source)`.

`_post_load_stack!(source, stack, nyears)` is called after loading so
sources can normalise sub-daily Ti before the canonical pipeline sees the data.
"""
function _load_weather(source::Type, area::Extent, years)
    primary_stack = _load_layers(weather_loader(source), source,
                                 primary_layers(source), area, years)
    baseline = fallback_source(source)
    stack = baseline === nothing ?
        RasterStack(primary_stack) :
        RasterStack(merge(primary_stack,
            _load_layers(weather_loader(baseline), baseline,
                         fallback_layers(source), area, years)))
    # Source files declare a `missingval`, so loaded eltypes are
    # `Union{Missing, T}`. The per-pixel reader does `value * unit`,
    # which throws `convert(Missing, Quantity)` on any masked cell.
    # `replace_missing(..., NaN)` strips the Missing union and lets
    # NaN propagate visibly if any cell is genuinely masked.
    stack = Rasters.replace_missing(stack, NaN)
    return _post_load_stack!(source, stack, length(years))
end

# ---------------------------------------------------------------------------
# Per-loader file-reading
# ---------------------------------------------------------------------------

function _load_layers(::YearlyTimeSeries, source, fields::Tuple, area::Extent, years)
    extras = _extra_getraster_kwargs(source)
    layers = map(fields) do name
        per_year = map(years) do yr
            path = getraster(source, name; date = Date(yr), extras...)
            read(crop(Raster(path; lazy = true); to = area, touches = true))
        end
        cat(per_year...; dims = Ti)
    end
    return NamedTuple{fields}(layers)
end

function _load_layers(::MonthlyClimatology, source, fields::Tuple, area::Extent, years)
    nyears = length(years)
    layers = map(fields) do name
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

function _load_layers(::DailyFiles, source, fields::Tuple, area::Extent, years)
    extras = _extra_getraster_kwargs(source)
    dates = _daily_date_sequence(years)
    nsteps = length(dates)
    layers = map(fields) do name
        per_day = map(dates) do d
            path = getraster(source, name; date = d, extras...)
            read(crop(Raster(path; lazy = true); to = area, touches = true))
        end
        first_day = first(per_day)
        spatial_dims = dims(first_day)
        stacked = cat(map(parent, per_day)...; dims = 3)
        Raster(stacked, (spatial_dims..., Ti(1:nsteps)); crs = crs(first_day))
    end
    return NamedTuple{fields}(layers)
end

# 365-day-per-year date sequence, dropping Feb 29 in leap years so the
# loaded data matches the fixed `nyears * 365` buffer.
function _daily_date_sequence(years)
    dates = Date[]
    for y in years
        for d in Date(y, 1, 1):Day(1):Date(y, 12, 31)
            (Dates.month(d) == 2 && Dates.day(d) == 29) && continue
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
    # `view(... X(lo..hi), Y(lo..hi), Ti(t1..t2))` selects by dim name —
    # storage order in the Zarr store (typically `(time, lat, lon)`) is
    # irrelevant since downstream access uses `[I..., Ti(k)]` dim wrappers.
    layers = map(fields) do name
        long_name = layername(source, name)
        raw = getproperty(full_stack, Symbol(long_name))
        read(view(raw,
            X(area.X[1] .. area.X[2]),
            Y(area.Y[1] .. area.Y[2]),
            Ti(time_start .. time_end),
        ))
    end
    return NamedTuple{fields}(layers)
end

# Build one 3-D climatology raster (X, Y, Ti = 12 * nyears) by reading 12
# monthly 2-D rasters via `path_for(month)`, stacking along a new Ti dim,
# then tiling across `nyears`. Shared by the climatology loaders.
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

# Unified derivation chain — every derivation any source might need, in
# topological order. Each step is skipped when its output is in the
# source's native variable set (so e.g. CHELSA{Climate} skips both the
# RH derivations and `cloud_cover`, WorldClim skips `actual_vapour_pressure`,
# TerraClimate skips both `wind_speed` and `vapour_pressure_deficit`).
#
# Derivations whose output is *not* native AND whose inputs are zero
# (because no upstream step populates them) compute harmless garbage that
# nothing downstream reads. We rely on the unused-output property rather
# than runtime input-availability checks.
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

# Hourly-mode chain — for sources like ERA5 that provide hourly values
# directly. No min/max envelope to build; the canonical hourly buffers
# *are* the env_hourly arrays. Just convert components into the shapes
# the solver wants and aggregate to env_daily.
const _DERIVATIONS_HOURLY = (
    Val(:wind_speed), # u/v → reference_wind_speed
    Val(:actual_vapour_pressure_from_dewpoint), # T_dew → actual_VP
    Val(:reference_humidity), # actual_VP + T → RH
    Val(:rainfall_daily_from_hourly), # hourly rainfall → daily total
    Val(:deep_soil_temperature_from_hourly), # annual mean of hourly T
)

# Sub-daily (6-hourly) chain — for sources like NCEPHourly that provide
# 6-hourly values that must be interpolated/disaggregated to hourly.
#
#  1. :solar_geometry          — populate scratch.solar.out with hourly clear-sky
#  2. :wind_speed_6h           — u/v → scalar wind at 6h resolution
#  3. :actual_vapour_pressure_6h — specific_humidity + pressure → VP at 6h
#  4. :interpolate_met_6h_to_1h  — linear 6h→1h for T, wind, VP, pressure, LW
#  5. :disaggregate_radiation_6h_to_1h — solar-aware 6h→1h for shortwave
#  6. :reference_humidity      — actual_VP + T → RH (reuses hourly-mode method)
#  7. :rainfall_daily_from_6h  — sum 4 × 6h blocks per day → daily total
#  8. :deep_soil_temperature_from_hourly — annual mean (reuses hourly-mode method)
const _DERIVATIONS_6H_TO_1H = (
    Val(:solar_geometry),
    Val(:wind_speed_6h),
    Val(:actual_vapour_pressure_6h),
    Val(:interpolate_met_6h_to_1h),
    Val(:disaggregate_radiation_6h_to_1h),
    Val(:reference_humidity),
    Val(:rainfall_daily_from_6h),
    Val(:deep_soil_temperature_from_hourly),
)

"""
    weather_derivations(::Type{<:RasterDataSource}) -> Tuple{Val, …}

The ordered derivation chain to run for this source. Defaults to
`_DEFAULT_DERIVATIONS`, which covers every monthly/daily source we
currently support via the skip-if-native logic. Hourly sources should
override to `_DERIVATIONS_HOURLY`.
"""
@inline weather_derivations(::Type) = _DEFAULT_DERIVATIONS

# ---------------------------------------------------------------------------
# Buffer allocation
# ---------------------------------------------------------------------------

"""
    allocate_weather_buffers(source, nyears) -> NamedTuple

Per-worker scratch for `assemble_weather!`. Allocates one unit-tagged
buffer per canonical variable plus the pre-constructed env structs that
share storage with those buffers.

The total number of timesteps is `nyears * steps_per_year(temporal_resolution(source))`,
and the env-minmax struct type is `MonthlyMinMaxEnvironment` or
`DailyMinMaxEnvironment` per the source's `temporal_resolution`. Everything
else — `DailyTimeseries`, `HourlyTimeseries`, the intermediate buffers —
just scales with the total timestep count.
"""
function allocate_weather_buffers(source::Type, days_of_year::AbstractVector{Int})
    _allocate_weather_buffers(temporal_resolution(source), source, days_of_year)
end

function _allocate_weather_buffers(resolution::Union{MonthlyResolution, DailyResolution},
                                   ::Any, days_of_year::AbstractVector{Int})
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

    # `_minmax_env_type(resolution)` returns the type constructor — same
    # field set for Monthly and Daily, so the kwargs are identical.
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
    # Pressure is constant per pixel — a `Fill` is rebuilt and swapped in
    # via `setproperties` each call, so this initial value is just a stub.
    # HourlyTimeseries is always 24× the main step count.
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
        days_of_year,
        environment_minmax, environment_daily, environment_hourly,
    )
end

# Hourly mode: env_minmax = nothing, env_hourly arrays ARE the canonical
# buffers (shared storage), env_daily arrays sit at 1/24 the resolution
# and are filled by aggregating the hourly canonical buffers.
function _allocate_weather_buffers(::HourlyResolution, ::Any, days_of_year::AbstractVector{Int})
    ndays  = length(days_of_year)    # one per unique day (Feb 29 dropped)
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

    # `environment_daily` lives at one entry per unique day — aggregated from hourly.
    rainfall_daily = zeros(typeof(0.0u"kg/m^2"), ndays)
    deep_soil_temperature = zeros(typeof(0.0u"K"), ndays)

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
        dewpoint_temperature, actual_vapour_pressure,
        deep_soil_temperature,
        days_of_year,
        environment_minmax = nothing,
        environment_daily, environment_hourly,
    )
end

# Sub-daily (6-hourly) mode: two tiers of buffers.
#   * `*_6h` staging buffers at native 6h resolution (1460 steps/year) — written by
#     `_read_native!` via the source's WeatherVariable declarations.
#   * Hourly output buffers (8760 steps/year) — filled by `_DERIVATIONS_6H_TO_1H`
#     and shared with `environment_hourly` exactly as in `HourlyResolution`.
#   * Daily aggregate buffers (365/year) shared with `environment_daily`.
function _allocate_weather_buffers(::SixHourlyResolution, _source,
                                   days_of_year::AbstractVector{Int})
    ndays  = length(days_of_year)   # one per unique day (Feb 29 dropped)
    n6h    = ndays * 4              # native 6h staging steps
    nhours = ndays * 24            # hourly output
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
        rainfall_daily, deep_soil_temperature,
        days_of_year,
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
)
    buffers = scratch.weather
    atm_pressure = atmospheric_pressure(site.elevation)
    variables = weather_variables(source)
    spy = steps_per_year(temporal_resolution(source))

    ctx = (;
        site, grid_elevation, lapse_rate_model, vapour_pressure_method,
        atmospheric_pressure = atm_pressure, scratch,
        steps_per_year = spy,
    )

    _read_native!(buffers, weather, variables, I)
    _run_derivations!(buffers, ctx, source, variables, canonical_overrides)
    _apply_canonical_overrides!(buffers, canonical_overrides, I)

    environment_hourly = _finalize_environment_hourly(buffers, variables, atm_pressure)

    return (; buffers.environment_minmax, buffers.environment_daily, environment_hourly)
end

# In monthly/daily mode, pressure is constant per pixel and the env_hourly
# struct's `pressure` field is a `Fill` derived from site elevation —
# rebuild and swap in. In hourly/sub-daily mode the pressure array is either
# loaded natively (ERA5 `:sp`) or interpolated from 6h staging buffers
# (NCEPHourly `surface_pressure_6h → pressure`), and must NOT be overwritten.
#
# The guard is `environment_minmax === nothing`: that flag is set exclusively
# by hourly/sub-daily allocators, so it unambiguously identifies modes where
# `environment_hourly.pressure` is already a populated concrete array.
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

# Override application — dispatch on Ti presence so both modes share one
# code path:
#   - constant-in-time override: spatial-only Raster, rank == length(I).
#     Grid mode: 2-D (X, Y); points mode: 1-D (Dim{:point},).
#   - time-varying override: trailing Ti dim, rank == length(I) + 1.
#     Grid mode: 3-D (X, Y, Ti); points mode: 2-D (Dim{:point}, Ti).
# `hasdim(raster, Ti)` is a compile-time check on the Raster's dim tuple,
# so the branch folds away for any concrete Raster type.
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
        @. out = lapse_adjust_temperature(lapse_rate_model, src, Δz)
    end
    return nothing
end

function derive!(::Val{:mean_temperature}, buffers, ctx)
    @. buffers.mean_temperature =
        (buffers.reference_temperature_max + buffers.reference_temperature_min) / 2
    return nothing
end

# Used by sources (NCEP) that provide specific humidity + surface pressure
# rather than vapour pressure directly:
#   actual_vapour_pressure = q · p / (0.622 + 0.378 · q)
# where q is specific humidity (kg/kg, dimensionless) and p is surface
# pressure. Result has units of pressure; auto-converted to kPa on
# assignment to the canonical buffer.
function derive!(::Val{:actual_vapour_pressure}, buffers, ctx)
    @inbounds for k in eachindex(buffers.actual_vapour_pressure)
        q = buffers.specific_humidity[k]
        p = buffers.surface_pressure[k]
        buffers.actual_vapour_pressure[k] = q * p / (0.622 + 0.378 * q)
    end
    return nothing
end

# Used by sources that provide actual vapour pressure (e.g. WorldClim's
# `:vapr` → `actual_vapour_pressure`, or NCEP via the SH→VP derivation
# above). Converts it into the VPD that the RH derivations consume:
#   VPD = e_s(mean_T) - actual_VP.
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

# 10 m → 2 m wind power-law (exponent 0.15).
const _WIND_10M_TO_2M_SHEAR = (2.0 / 10.0)^0.15

function derive!(::Val{:reference_wind_max}, buffers, ctx)
    @. buffers.reference_wind_max = buffers.wind_speed * _WIND_10M_TO_2M_SHEAR
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
    spy = ctx.steps_per_year
    nyears = length(buffers.mean_temperature) ÷ spy
    _annual_means!(buffers.deep_soil_temperature,
                   buffers.mean_temperature, nyears, spy)
    return nothing
end

# ---------------------------------------------------------------------------
# Hourly-mode derivations
# ---------------------------------------------------------------------------

# Actual vapour pressure from 2 m dewpoint temperature (ERA5):
# at the dewpoint, the air is saturated → actual vapour pressure equals
# saturation vapour pressure at the dewpoint.
function derive!(::Val{:actual_vapour_pressure_from_dewpoint}, buffers, ctx)
    method = ctx.vapour_pressure_method
    @inbounds for k in eachindex(buffers.actual_vapour_pressure)
        buffers.actual_vapour_pressure[k] =
            vapour_pressure(method, buffers.dewpoint_temperature[k])
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

# Aggregate 24 hourly rainfall values into one daily total (kg/m²) for
# env_daily.rainfall.
# TODO: propagate hourly rainfall into Microclimate.jl. The solver should
# consume the hourly array (`env_hourly.rainfall`) directly at the bottom
# level rather than us pre-aggregating to daily totals here — that loses
# sub-daily intensity information that affects soil infiltration / runoff.
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

# Annual mean of hourly reference_temperature, broadcast across all 365
# daily entries of `deep_soil_temperature`.
function derive!(::Val{:deep_soil_temperature_from_hourly}, buffers, ctx)
    ndays  = length(buffers.deep_soil_temperature)
    nhours = length(buffers.reference_temperature)
    if ndays < 365
        # Sub-annual run: fill with mean of all available hours.
        run_mean = sum(buffers.reference_temperature) / nhours
        fill!(buffers.deep_soil_temperature, run_mean)
        return nothing
    end
    nyears = ndays ÷ 365
    @inbounds for y in 1:nyears
        h0 = (y - 1) * 8760
        annual_sum = zero(eltype(buffers.reference_temperature))
        for h in 1:8760
            annual_sum += buffers.reference_temperature[h0 + h]
        end
        annual_mean = annual_sum / 8760
        d0 = (y - 1) * 365
        for d in 1:365
            buffers.deep_soil_temperature[d0 + d] = annual_mean
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Sub-daily (6-hourly) derivations
# ---------------------------------------------------------------------------

# Run solar_radiation! to populate scratch.solar.out with hourly clear-sky
# global_horizontal values for all days. Must be the first step in
# _DERIVATIONS_6H_TO_1H so subsequent steps can read the clear-sky output.
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

# Scalar wind speed from U/V components at 6h resolution.
function derive!(::Val{:wind_speed_6h}, buffers, _ctx)
    @. buffers.wind_speed_6h = sqrt(buffers.u_wind_6h^2 + buffers.v_wind_6h^2)
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

# Linear interpolation helper: 1460 (4×365) 6-hourly values → 8760 hourly.
# The four 6h blocks per day have their midpoints at hours 3, 9, 15, 21.
# Values outside the first/last block boundary are clamped to the edge value.
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

# Linearly interpolate met variables (T, wind, humidity, pressure, longwave)
# from 6h to hourly. Uses _interp_6h_to_1h! for each variable.
function derive!(::Val{:interpolate_met_6h_to_1h}, buffers, _ctx)
    _interp_6h_to_1h!(buffers.reference_temperature,    buffers.reference_temperature_6h)
    _interp_6h_to_1h!(buffers.wind_speed,               buffers.wind_speed_6h)
    _interp_6h_to_1h!(buffers.actual_vapour_pressure,   buffers.actual_vapour_pressure_6h)
    _interp_6h_to_1h!(buffers.pressure,                 buffers.surface_pressure_6h)
    _interp_6h_to_1h!(buffers.longwave_radiation,       buffers.longwave_radiation_6h)
    return nothing
end

# Solar-geometry-aware shortwave disaggregation from 6h to hourly.
#
# Algorithm (NicheMapR/microclima approach):
#   For each 6h block, compute the clear-sky mean from scratch.solar.out.global_horizontal
#   (populated by derive!(:solar_geometry)), then derive opacity = observed / clear_sky.
#   Assign constant opacity within the block and multiply by the hourly clear-sky.
#   Days where the clear-sky block mean is zero (polar night) get zero radiation.
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

Derive monthly relative humidity following NicheMapR's micro_terra.R:
  actual_vapour_pressure = saturation_vapour_pressure(mean_temperature)
                           − vapour_pressure_deficit
  relative_humidity = actual_vapour_pressure
                           / saturation_vapour_pressure(reference_temperature)

Pass `reference_temperature = minimum_temperature` for RH_max and
`reference_temperature = maximum_temperature` for RH_min. Output is clamped
to [0, 1].
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
    _annual_means!(out, values, nyears, steps_per_year)

For each of `nyears` years (`steps_per_year` consecutive entries), every
entry within the year gets replaced by that year's annual mean. Used to
build `deep_soil_temperature` for both monthly (`steps_per_year = 12`) and
daily (`steps_per_year = 365`) resolutions.
"""
function _annual_means!(out::AbstractVector,
                        values::AbstractVector,
                        nyears::Int,
                        steps_per_year::Int)
    if nyears == 0
        # Sub-annual run: fill with mean of all available steps.
        run_mean = sum(values) / length(values)
        fill!(out, run_mean)
        return out
    end
    @inbounds for y in 1:nyears
        annual_sum = zero(eltype(values))
        for k in 1:steps_per_year
            annual_sum += values[(y - 1) * steps_per_year + k]
        end
        annual_mean = annual_sum / steps_per_year
        for k in 1:steps_per_year
            out[(y - 1) * steps_per_year + k] = annual_mean
        end
    end
    return out
end
