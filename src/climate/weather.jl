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
the source natively provides. A physics derivation runs only if its output
is not already in this native set (`_is_native`).
"""
function weather_variables end

# ---------------------------------------------------------------------------
# Time: two orthogonal dials — calendar and timestep
# ---------------------------------------------------------------------------
#
# A run's time is two independent choices, never one ladder:
#
#   * calendar — which days we run. A month is either one representative
#     day (`Monthly`) or all of its days (`Daily`). This is the long
#     accounting axis.
#   * timestep — the step *within* a day. Either a per-step min/max envelope
#     the solver expands into a diel curve (`MinMax`), or `N` evenly spaced
#     samples per day (`SubDaily{N}` — hourly is `SubDaily{24}`, six-hourly
#     `SubDaily{4}`). This is the sub-daily axis.
#
# They do not trade off: hourly runs inside each day whether the month
# contributes one day or thirty. A source declares the calendar it covers
# (`weather_calendar`) and the timestep its data arrives at (`native_timestep`).
# The run declares the timestep it wants out (the target, default hourly).
# Moving a quantity between timesteps is `resample!`, driven by each
# quantity's own `resampling_rule` — never by which source it came from.

abstract type Calendar end

"""
    Monthly

12 representative mid-month days per simulated year. Used by TerraClimate,
CHELSA, WorldClim. Picks `MonthlyMinMaxEnvironment`, `Microclimate.DEFAULT_DAYS`
for solar geometry, and non-consecutive-day mode (each month an independent
representative day).
"""
struct Monthly <: Calendar end

"""
    Daily

Every day, 365 per simulated year (Feb 29 dropped). Used by daily and
sub-daily sources (NCEP, AWAP, ERA5, GRIDMET). Picks `DailyMinMaxEnvironment`
(in `MinMax` timestep), `1:365` for solar geometry, and consecutive-day mode
(each day inherits soil state from the previous).
"""
struct Daily <: Calendar end

abstract type Timestep end

"""
    MinMax

The source provides a per-step min/max envelope and no within-day samples;
the solver synthesises the diel curve. Used by monthly/daily sources.
"""
struct MinMax <: Timestep end

"""
    SubDaily{N}

`N` evenly spaced samples within each day. `Hourly == SubDaily{24}`,
`SixHourly == SubDaily{4}`. A quantity moves between sub-daily timesteps (and
to/from a day total) generically via `resample!` and its `resampling_rule`.
"""
struct SubDaily{N} <: Timestep end
const Hourly = SubDaily{24}
const SixHourly = SubDaily{4}

@inline samples_per_day(::SubDaily{N}) where {N} = N
@inline samples_per_day(::MinMax) = 1

"""
    weather_calendar(::Type{<:RasterDataSource}) -> Calendar

Which days the source covers. Default `Monthly()`; daily/sub-daily sources
override to `Daily()`.
"""
@inline weather_calendar(::Type) = Monthly()

"""
    native_timestep(::Type{<:RasterDataSource}) -> Timestep

The within-day timestep the source's data arrives at. Default `MinMax()`
(daily/monthly min/max sources); hourly/sub-daily sources override.
"""
@inline native_timestep(::Type) = MinMax()

# Representative (solar-geometry) days per year, by calendar.
@inline days_per_year(::Monthly) = 12
@inline days_per_year(::Daily)   = 365

# Total per-year steps at a given calendar × timestep. MinMax → one step per
# day (12 or 365); SubDaily{N} → N per day. Native timestep sizes the native
# buffers; target timestep sizes the output buffers.
@inline steps_per_year(cal::Calendar, cad::Timestep) =
    days_per_year(cal) * samples_per_day(cad)

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
    MonthlyClimatologyPeriod

12 files per layer for a climatology covering a specific period (a finite
window, e.g. a decade) rather than a fixed baseline. Like `MonthlyClimatology`
but with an extra `date` kwarg selecting that period's window. Used by
CHELSA{Future{Climate}}. `getraster(source, name; date = period_date, month)`.
"""
struct MonthlyClimatologyPeriod <: WeatherLoader end

"""
    MultiBandClimatologyPeriod

One multi-band file per layer for a specific period's date — the 12
months live inside a single GeoTIFF as separate bands. Used by
WorldClim{Future{Climate}}. `getraster(source, name; date = period_date)`.
"""
struct MultiBandClimatologyPeriod <: WeatherLoader end

"""
    DailyFiles

One file per day per layer — the source ships one 2-D raster per calendar
day. Used by daily station/grid products such as AWAP.
`getraster(source, name; date = day)`. Feb 29 in leap years is skipped
to match the fixed `nyears * 365` buffer size.
"""
struct DailyFiles <: WeatherLoader end

"""
    ContiguousTimeSeries

The whole time series in a single store — every layer with one unbroken `Ti`
axis — rather than partitioned into per-period files. The run's X/Y/time
window is read by slicing that store; no per-date iteration. Used by ERA5,
whose `getraster(source)` (no args) returns a `CachedCloudSource` opened once
through Rasters, each layer read by its long name
(`RasterDataSources.layername(source, sym)`). The storage backend (Zarr, via
the ARCO-ERA5 store today) is irrelevant — Rasters abstracts it, and any
format organised this way loads identically.
"""
struct ContiguousTimeSeries <: WeatherLoader end

"""
    weather_loader(::Type{<:RasterDataSource}) -> WeatherLoader

Singleton declaring how this source's files are organised. One of
`YearlyTimeSeries`, `MonthlyClimatology`, `MonthlyClimatologyPeriod`,
`MultiBandClimatologyPeriod`.
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
call (in addition to the loader's own kwargs like `date`/`month`). Default empty.
"""
@inline _extra_getraster_kwargs(::Type) = (;)

"""
    weather_grid_elevation(source, weather, I) -> Quantity or nothing

Return the weather grid's elevation at spatial index `I` as a `Unitful` length
(e.g. `450.0u"m"`), or `nothing` if the source does not carry its own elevation.

When `nothing` is returned, lapse rate correction is skipped (effectively
`grid_elevation = site.elevation`, so `Δz = 0`). Sources that ship their own
elevation layer (e.g. CRUCL2 via its `elv` variable) should override this to
activate the lapse correction machinery in `_lapse_correct!`.
"""
@inline weather_grid_elevation(::Type, _weather, _I) = nothing

# ---------------------------------------------------------------------------
# Generic _load_weather
# ---------------------------------------------------------------------------

"""
    _post_load_stack!(source, stack, nyears) -> RasterStack

Source-specific post-processing hook called immediately after every layer of
`stack` has been loaded. Default is a no-op (returns `stack` unchanged).

Override for sources whose files contain a finer Ti than the native buffers
expect and need aggregating after load (see `_aggregate_ti_to_daily`).
"""
_post_load_stack!(::Type, stack, _nyears) = stack

# Average a sub-daily Ti layer down by `factor` steps using Rasters.aggregate.
function _aggregate_ti_to_daily(layer::AbstractRaster, factor::Int)
    return Rasters.aggregate(mean, layer, (Ti(factor),))
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
    stack = _post_load_stack!(source, stack, length(years))
    # Every source must deliver exactly the native step count it declares
    # (calendar days/year × samples/day × years). Catch any loader that
    # returns a different Ti rather than mis-slicing it downstream.
    expected = steps_per_year(weather_calendar(source), native_timestep(source)) * length(years)
    n = length(dims(stack, Ti))
    n == expected || error(
        "$(nameof(source)): loaded Ti=$n, but the source declares $expected " *
        "($(nameof(typeof(weather_calendar(source)))) × " *
        "$(nameof(typeof(native_timestep(source)))), $(length(years)) years)")
    return stack
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

function _load_layers(::MonthlyClimatologyPeriod, source, fields::Tuple, area::Extent, years)
    nyears = length(years)
    period_date = Date(first(years))
    layers = map(fields) do name
        _build_monthly_climatology(area, nyears) do month
            getraster(source, name; date = period_date, month)
        end
    end
    return NamedTuple{fields}(layers)
end

function _load_layers(::MultiBandClimatologyPeriod, source, fields::Tuple,
                      area::Extent, years)
    nyears = length(years)
    period_date = Date(first(years))
    layers = map(fields) do name
        path = getraster(source, name; date = period_date)
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

function _load_layers(::ContiguousTimeSeries, source, fields::Tuple, area::Extent, years)
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

# Physics chains — ordered lists of timestep-blind derivations. A chain is a
# property of the solver-input shape it feeds, not of any source: `_is_native`
# skips whatever the source already provides (so CHELSA{Climate} skips the
# relative-humidity derivations and `cloud_cover`, TerraClimate skips
# `wind_speed`, etc.). Time
# never appears here — moving a quantity between timesteps is `resample!`.
#
# A derivation whose output is not native AND whose inputs are unpopulated
# computes harmless garbage nothing downstream reads; we rely on that rather
# than runtime input-availability checks.

# Envelope shape (daily/monthly min-max sources): builds the min/max
# environment the solver expands into a diel curve.
const _ENVELOPE_PHYSICS = (
    Val(:reference_temperature_max), # ← lapse(maximum_temperature)
    Val(:reference_temperature_min), # ← lapse(minimum_temperature)
    Val(:mean_temperature),          # ← (ref_max + ref_min) / 2
    Val(:wind_speed),                # ← √(eastward_wind² + northward_wind²)
    Val(:actual_vapour_pressure),    # ← q·p / (0.622 + 0.378·q)
    Val(:vapour_pressure_deficit),   # ← saturation_vapour_pressure(mean_temperature) − actual_vapour_pressure
    Val(:reference_humidity_max),    # ← from vapour_pressure_deficit + ref_temp_min
    Val(:reference_humidity_min),    # ← from vapour_pressure_deficit + ref_temp_max
    Val(:reference_wind_max),        # ← wind_speed × shear_factor
    Val(:reference_wind_min),        # ← ref_wind_max × 0.1
    Val(:cloud_cover),               # ← from global_radiation
    Val(:cloud_min),                 # ← cloud_cover × 0.5 clamp
    Val(:cloud_max),                 # ← cloud_cover × 2.0 clamp
    Val(:deep_soil_temperature),     # ← annual mean of mean_temperature
)

# Series shapes (within-day samples). The native combiners run on the native
# tier (before resampling); the output physics run on the output tier (after).
# `solar_geometry` prepares the clear-sky baseline the shortwave resample reads.
const _NATIVE_PHYSICS = (
    Val(:solar_geometry),         # populate clear-sky for the shortwave rule
    Val(:wind_speed),             # u/v → wind speed
    Val(:actual_vapour_pressure), # specific-humidity·pressure or dewpoint → actual vapour pressure
)
const _OUTPUT_PHYSICS = (
    Val(:apply_wind_shear),       # measurement height → reference height
    Val(:reference_humidity),     # actual_vapour_pressure + temperature → relative humidity
    Val(:deep_soil_temperature),  # annual mean of air temperature
    Val(:rainfall_daily),         # accumulate the finest rainfall → daily total
)

# ---------------------------------------------------------------------------
# Temporal resampling — the time axis, generic over timestep
# ---------------------------------------------------------------------------
#
# `resample!` moves a quantity from a native timestep (`nin` samples/day) to a
# target (`nout` samples/day, or a daily total). How a quantity resamples is a
# property of the quantity, declared once by `resampling_rule`:
#   * Interpolate  — densify a smooth field (temperature, wind, vapour
#                    pressure, pressure, longwave).
#   * Disaggregate — split a block aggregate onto a reference profile: hold
#                    opacity constant within each native block and multiply by
#                    the clear-sky curve (shortwave).
#   * Accumulate   — fluxes that sum over a day (rainfall).

abstract type ResamplingRule end
struct Interpolate <: ResamplingRule end
struct Disaggregate <: ResamplingRule end
struct Accumulate <: ResamplingRule end

@inline resampling_rule(::Val) = Interpolate()
@inline resampling_rule(::Val{:global_radiation}) = Disaggregate()
@inline resampling_rule(::Val{:rainfall}) = Accumulate()

# Native sample `i` (0-based) of day `d`, clamping across day boundaries to the
# neighbouring day's edge sample (or this day's edge at the series ends).
@inline function _native_at(src, d, i, nin, ndays)
    if i < 0
        d > 1 ? src[(d - 2) * nin + nin] : src[(d - 1) * nin + 1]
    elseif i >= nin
        d < ndays ? src[d * nin + 1] : src[(d - 1) * nin + nin]
    else
        src[(d - 1) * nin + i + 1]
    end
end

# Linear interpolation from `nin` to `nout` samples/day. Each sample sits at the
# midpoint of its slot; output midpoints interpolate between bracketing native
# midpoints. Reduces to the old 6h→1h kernel for nin=4, nout=24.
function _resample!(out, src, ::Interpolate, nin::Int, nout::Int)
    ndays = length(src) ÷ nin
    @inbounds for d in 1:ndays
        ob = (d - 1) * nout
        for j in 0:(nout - 1)
            x  = ((j + 0.5) / nout) * nin - 0.5   # output midpoint on native index axis
            i0 = floor(Int, x)
            t  = x - i0
            v0 = _native_at(src, d, i0,     nin, ndays)
            v1 = _native_at(src, d, i0 + 1, nin, ndays)
            out[ob + j + 1] = v0 * (1 - t) + v1 * t
        end
    end
    return out
end

# Solar-geometry-aware shortwave disaggregation: per native block, opacity =
# observed / clear-sky block-mean, applied to each output step's clear-sky.
# `clear_sky` is the hourly (nout/day) clear-sky global-horizontal series.
function _resample!(out, src, ::Disaggregate, nin::Int, nout::Int, clear_sky)
    ndays = length(src) ÷ nin
    hpb = nout ÷ nin                       # output steps per native block
    @inbounds for d in 1:ndays
        hb = (d - 1) * nout
        for b in 0:(nin - 1)
            cs_sum = zero(eltype(clear_sky))
            for j in 1:hpb
                cs_sum += clear_sky[hb + b * hpb + j]
            end
            cs_mean = cs_sum / hpb
            obs = src[(d - 1) * nin + b + 1]
            opacity = cs_mean <= zero(cs_mean) ? 0.0 :
                clamp(ustrip(u"W/m^2", obs) / ustrip(u"W/m^2", cs_mean), 0.0, 1.0)
            for j in 1:hpb
                out[hb + b * hpb + j] = opacity * clear_sky[hb + b * hpb + j]
            end
        end
    end
    return out
end

# Sum the `nin` native samples in each day into one daily total.
function _resample!(daily, src, ::Accumulate, nin::Int)
    ndays = length(daily)
    @inbounds for d in 1:ndays
        s = zero(eltype(src))
        for i in 1:nin
            s += src[(d - 1) * nin + i]
        end
        daily[d] = s
    end
    return daily
end

# Carry every output quantity that also exists in the native tier from native to
# output, each by its own `resampling_rule`. Output-only quantities (no native
# source) and `Accumulate` quantities (which target the daily tier, via
# `derive!(:rainfall_daily)`) are left to the physics chains. Driven entirely by
# the rule trait — no quantity is named here.
@inline function _resample_to_output!(buffers, nin::Int, nout::Int, ctx)
    native = buffers.native
    map(_SERIES_OUTPUT) do name
        _resample_quantity!(buffers, native, Val(name), nin, nout, ctx)
    end
    return nothing
end

@inline function _resample_quantity!(buffers, native, ::Val{name}, nin, nout, ctx) where {name}
    hasproperty(native, name) || return nothing
    _resample_by_rule!(getproperty(buffers, name), getproperty(native, name),
                       resampling_rule(Val(name)), nin, nout, ctx)
    return nothing
end

# Interpolate and Disaggregate reconstruct sub-steps the native tier
# doesn't carry; reconstruction is a coarse→fine operation. At equal cadence
# there is nothing to reconstruct, so the carry is identity — for every
# quantity. This is why a source whose native timestep already matches the
# target (ERA5 hourly → hourly) keeps its observed values untouched instead of
# being routed through a disaggregation kernel built for the coarse case.
@inline function _resample_by_rule!(out, src, rule::Union{Interpolate,Disaggregate},
                                    nin, nout, ctx)
    nin == nout ? copyto!(out, src) : _disaggregate!(out, src, rule, nin, nout, ctx)
    return nothing
end
@inline _disaggregate!(out, src, rule::Interpolate, nin, nout, ctx) =
    _resample!(out, src, rule, nin, nout)
@inline _disaggregate!(out, src, rule::Disaggregate, nin, nout, ctx) =
    _resample!(out, src, rule, nin, nout, ctx.scratch.solar.out.global_horizontal)
# Accumulate quantities target the daily tier, not the output tier — handled by
# `derive!(:rainfall_daily)`, so skip them here.
@inline _resample_by_rule!(out, src, ::Accumulate, nin, nout, ctx) = nothing

# ---------------------------------------------------------------------------
# Shared physics combiners — one implementation, called at any timestep
# ---------------------------------------------------------------------------

# Scalar wind speed from horizontal U/V components.
@inline function _wind_speed!(out, u, v)
    @. out = sqrt(u^2 + v^2)
    return out
end

# Actual vapour pressure from specific humidity (kg/kg) and pressure:
#   e = q·p / (0.622 + 0.378·q)
function _vapour_pressure_from_specific_humidity!(out, q, p)
    @inbounds for k in eachindex(out)
        out[k] = q[k] * p[k] / (0.622 + 0.378 * q[k])
    end
    return out
end

# Annual mean of a temperature `series`, replicated across each year's entries of
# `out`. `out_per_year` is the number of `out` entries per simulated year (12
# monthly, 365 daily); the year count comes from that, and the series is split
# into the same number of years. A sub-annual run collapses to one block (mean of
# everything). Serves both the envelope mean (out and series same length) and the
# hourly series (series is 24× out).
function _deep_soil!(out, series, out_per_year)
    nyears = max(1, length(out) ÷ out_per_year)
    out_per = length(out) ÷ nyears
    series_per = length(series) ÷ nyears
    @inbounds for y in 1:nyears
        s = zero(eltype(series))
        sb = (y - 1) * series_per
        for k in 1:series_per
            s += series[sb + k]
        end
        m = s / series_per
        ob = (y - 1) * out_per
        for d in 1:out_per
            out[ob + d] = m
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Buffer allocation
# ---------------------------------------------------------------------------

"""
    allocate_weather_buffers(source, target_timestep, nyears) -> NamedTuple

Per-worker scratch for `assemble_weather!`. Allocates one unit-tagged
buffer per canonical variable plus the pre-constructed env structs that
share storage with those buffers.

Buffer sizes follow the two dials: `weather_calendar(source)` sets the days
per year (and, in `MinMax` timestep, the env-minmax struct —
`MonthlyMinMaxEnvironment` vs `DailyMinMaxEnvironment`); `native_timestep(source)`
sizes the native tier and `target_timestep` the output tier.
"""
function allocate_weather_buffers(source::Type, target::Timestep, nyears::Int)
    cal = weather_calendar(source)
    native = native_timestep(source)
    _allocate_weather_buffers(cal, native, target, _native_buffer_names(source),
                              _days_of_year(cal, nyears))
end

# The native tier holds every quantity read from the stack plus the two combiners
# derived at the native timestep before resampling. Built from the source's own
# `weather_variables`, so any sub-daily source rides the one path — there is no
# fixed input list and no separate "native == target" shape.
@inline _native_buffer_names(source::Type) =
    (map(canonical_name, weather_variables(source))..., :wind_speed, :actual_vapour_pressure)

# The working unit each canonical quantity is computed in — the single source
# of truth, so the allocators below just list quantity names. `nothing` (the
# default) means a dimensionless Float64 quantity (humidity fractions, cloud
# fractions, specific humidity, soil moisture). Distinct from `canonical_unit`,
# which is the unit the *output* layers are presented in (e.g. K here vs °C out).
working_units(::Val{:maximum_temperature})          = u"K"
working_units(::Val{:minimum_temperature})          = u"K"
working_units(::Val{:mean_temperature})             = u"K"
working_units(::Val{:reference_temperature})        = u"K"
working_units(::Val{:reference_temperature_min})    = u"K"
working_units(::Val{:reference_temperature_max})    = u"K"
working_units(::Val{:dewpoint_temperature})         = u"K"
working_units(::Val{:deep_soil_temperature})        = u"K"
working_units(::Val{:eastward_wind})                = u"m/s"
working_units(::Val{:northward_wind})               = u"m/s"
working_units(::Val{:wind_speed})                   = u"m/s"
working_units(::Val{:reference_wind_min})           = u"m/s"
working_units(::Val{:reference_wind_max})           = u"m/s"
working_units(::Val{:pressure})                     = u"Pa"
working_units(::Val{:actual_vapour_pressure})       = u"kPa"
working_units(::Val{:vapour_pressure_deficit})      = u"kPa"
working_units(::Val{:global_radiation})             = u"W/m^2"
working_units(::Val{:longwave_radiation})           = u"W/m^2"
working_units(::Val{:rainfall})                     = u"kg/m^2"
working_units(::Val{:rainfall_daily})               = u"kg/m^2"
working_units(::Val{:zenith_angle})                 = u"°"

@inline _zeros(unit, n) = zeros(typeof(0.0 * unit), n)
@inline _zeros(::Nothing, n) = zeros(Float64, n)
@inline _zero_buffer(::Val{name}, n) where {name} = _zeros(working_units(Val(name)), n)

# A NamedTuple of zeroed canonical buffers of length `n`, one per name, each
# tagged with `working_units`. (`names` is a literal tuple at every call site, so
# the per-name unit dispatch folds to a concrete buffer type.)
@inline _zero_buffers(names::Tuple, n) =
    NamedTuple{names}(map(name -> _zero_buffer(Val(name), n), names))

# `environment_daily` is the same per-step struct in every mode — only its
# length and the (shared) rainfall / deep-soil arrays differ.
function _environment_daily(n, rainfall, deep_soil_temperature)
    DailyTimeseries(;
        shade = zeros(Float64, n),
        soil_wetness = zeros(Float64, n),
        surface_emissivity = fill(0.95, n),
        cloud_emissivity = fill(0.95, n),
        rainfall, deep_soil_temperature,
        leaf_area_index = fill(0.1, n),
    )
end

# `environment_hourly` for the series modes: every field populated from the
# output-tier buffers `out`, which the solver reads via the `env_minmax::Nothing`
# dispatch.
_environment_hourly(out) = HourlyTimeseries(;
    out.pressure, out.reference_temperature, out.reference_humidity,
    reference_wind_speed = out.wind_speed,
    out.global_radiation, out.longwave_radiation, out.cloud_cover,
    out.rainfall, out.zenith_angle,
)

# `environment_hourly` for the envelope modes: a pressure-only stub — the solver
# synthesises the diel curve from env_minmax. Pressure is rebuilt/swapped per
# call, so the Fill value here is just a placeholder (24× the per-step count).
_environment_hourly_stub(n) = HourlyTimeseries(;
    pressure = Fill(atmospheric_pressure(0.0u"m"), n * 24),
    reference_temperature = nothing, reference_humidity = nothing,
    reference_wind_speed = nothing, global_radiation = nothing,
    longwave_radiation = nothing, cloud_cover = nothing,
    rainfall = nothing, zenith_angle = nothing,
)

# The per-step buffer set every envelope (min/max) source allocates — a superset:
# a source provides some quantities natively and derives the rest; any it never
# touches stay zero.
const _ENVELOPE_BUFFERS = (
    :maximum_temperature, :minimum_temperature, :mean_temperature,
    :reference_temperature_min, :reference_temperature_max,
    :wind_speed, :eastward_wind, :northward_wind, :reference_wind_min, :reference_wind_max,
    :vapour_pressure_deficit, :actual_vapour_pressure,
    :specific_humidity, :pressure,
    :reference_humidity_min, :reference_humidity_max,
    :global_radiation, :cloud_cover, :cloud_min, :cloud_max,
    :rainfall, :deep_soil_temperature, :soil_moisture,
)

# Every series-shape source (hourly or resampled-to-hourly) shares the same
# output tier and the same daily-aggregate tier — only whether it carries a
# separate coarser native tier differs (resample vs not).
const _SERIES_OUTPUT = (
    :reference_temperature, :reference_humidity, :wind_speed, :pressure,
    :cloud_cover, :global_radiation, :longwave_radiation, :rainfall,
    :zenith_angle, :actual_vapour_pressure,
)
const _SERIES_DAILY = (:rainfall_daily, :deep_soil_temperature, :soil_moisture)

# MinMax timestep (daily/monthly high-low sources): the solver draws the diel
# curve from the per-step envelope, so the target timestep does not change the
# buffers. The calendar picks the env-minmax struct.
function _allocate_weather_buffers(calendar::Calendar, ::MinMax, ::Timestep,
                                   native_names, days_of_year::AbstractVector{Int})
    nsteps = length(days_of_year)
    # Some of these arrays are also held inside the env structs below; writes
    # via the canonical name show up in the struct (same array).
    b = _zero_buffers(_ENVELOPE_BUFFERS, nsteps)

    # `_minmax_env_type(calendar)` returns the type constructor — same field set
    # for Monthly and Daily, so the kwargs are identical.
    environment_minmax = _minmax_env_type(calendar)(;
        b.reference_temperature_min, b.reference_temperature_max,
        b.reference_wind_min, b.reference_wind_max,
        b.reference_humidity_min, b.reference_humidity_max,
        b.cloud_min, b.cloud_max,
        minima_times = (temp = 0, wind = 0, humidity = 1, cloud = 1),
        maxima_times = (temp = 1, wind = 1, humidity = 0, cloud = 0),
    )
    environment_daily = _environment_daily(nsteps, b.rainfall, b.deep_soil_temperature)
    environment_hourly = _environment_hourly_stub(nsteps)

    return (; b..., days_of_year,
        environment_minmax, environment_daily, environment_hourly)
end

# Every sub-daily source (whatever its native timestep) rides one shape: an output
# tier at the target timestep, a daily-aggregate tier, and a separate `native` tier
# (nested — same quantity names at the native length) holding the source's native
# inputs plus the combiners derived before resampling. `assemble_weather!` derives
# physics on the native tier, resamples each quantity to the output tier, then runs
# the output physics. When the native timestep equals the target the resample is an
# identity — a value of the resample, not a separate code path.
function _allocate_weather_buffers(::Calendar, native_step::SubDaily, target_step::SubDaily,
                                   native_names, days_of_year::AbstractVector{Int})
    ndays = length(days_of_year)                  # 365 per year (no-leap calendar)
    nnat  = ndays * samples_per_day(native_step)  # native steps
    nout  = ndays * samples_per_day(target_step)  # output steps

    output = _zero_buffers(_SERIES_OUTPUT, nout)
    daily  = _zero_buffers(_SERIES_DAILY, ndays)
    native = _zero_buffers(native_names, nnat)

    environment_daily  = _environment_daily(ndays, daily.rainfall_daily, daily.deep_soil_temperature)
    environment_hourly = _environment_hourly(output)

    return (; output..., native, daily..., days_of_year,
        environment_minmax = nothing, environment_daily, environment_hourly)
end

# ---------------------------------------------------------------------------
# Generic assembler
# ---------------------------------------------------------------------------

"""
    assemble_weather!(scratch, weather, source, site, I; kwargs...) -> NamedTuple

Build `(environment_minmax, environment_daily, environment_hourly)` for one
pixel of `weather` (a `RasterStack` loaded by `_load_weather(source, …)`).
Works in place into the scratch buffers held by `scratch.weather`.

Dispatches on `(native_timestep, target_timestep)` to one of three shapes
(`_assemble!`): envelope (min/max), hourly-series in place, or native tier
+ `resample!` to the target. Each reads the native variables, runs the
timestep-blind physics, and lays out the daily/hourly tiers the solver reads;
then a fresh pressure `Fill` is swapped in for envelope-mode sources.
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
    target_timestep::Timestep = Hourly(),
)
    buffers = scratch.weather
    atm_pressure = atmospheric_pressure(site.elevation)
    variables = weather_variables(source)
    cal = weather_calendar(source)
    native = native_timestep(source)

    ctx = (;
        site, grid_elevation, lapse_rate_model, vapour_pressure_method,
        atmospheric_pressure = atm_pressure, scratch,
        days_of_year = buffers.days_of_year,
        # Output day count per year, for the deep-soil annual mean.
        steps_per_year = days_per_year(cal),
        native_timestep = native, target_timestep,
        wind_reference_height,
    )

    _assemble!(native, target_timestep, buffers, weather, source, ctx, variables,
               canonical_overrides, I)
    _apply_canonical_overrides!(buffers, canonical_overrides, I)

    environment_hourly = _finalize_environment_hourly(buffers, variables, atm_pressure)

    return (; buffers.environment_minmax, buffers.environment_daily, environment_hourly)
end

# Three assembly shapes, dispatched on (native timestep, target timestep). Each
# reads native variables, runs physics (timestep-blind), and lays out the daily
# / hourly tiers the solver reads. Kind (which native variables, which physics)
# is orthogonal to time (the timesteps here); time-shuffling is `resample!`.

# Envelope sources (daily/monthly high-low): the solver draws the diel curve
# from the per-step min/max, so there is nothing to resample.
@inline function _assemble!(::MinMax, ::Timestep, buffers, weather, source, ctx,
                            variables, overrides, I)
    _read_native!(buffers, weather, variables, I)
    _run_physics!(buffers, ctx, _ENVELOPE_PHYSICS, variables, overrides)
    return nothing
end

# Every sub-daily source: read into the native tier, run the native combiners
# there, resample every shared quantity to the output tier by its rule, then the
# output physics. When the native timestep equals the target (e.g. ERA5 hourly),
# the resample is an identity — no special path.
function _assemble!(native_step::SubDaily, target_step::SubDaily, buffers, weather, source,
                    ctx, variables, overrides, I)
    _read_native!(buffers.native, weather, variables, I)
    _run_physics!(buffers.native, ctx, _NATIVE_PHYSICS, variables, overrides)
    _resample_to_output!(buffers, samples_per_day(native_step), samples_per_day(target_step), ctx)
    _run_physics!(buffers, ctx, _OUTPUT_PHYSICS, variables, overrides)
    return nothing
end

# In envelope (MinMax) mode, pressure is constant per pixel and the env_hourly
# struct's `pressure` field is a `Fill` derived from site elevation — rebuild
# and swap in. In sub-daily mode the pressure array is either loaded natively
# (ERA5 `:sp`) or resampled from the native tier (NCEP surface pressure), and
# must NOT be overwritten.
#
# The guard is `environment_minmax === nothing`: that flag is set exclusively
# by the sub-daily allocators, so it unambiguously identifies modes where
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

# Run a fixed physics chain, skipping any output the source provides natively
# or the user overrides. The chain is a property of the solver-input shape
# (envelope vs hourly series), not of the source — kind enters only through
# the skip-if-native check.
@inline function _run_physics!(buffers, ctx, chain::Tuple, variables::Tuple,
                               overrides::NamedTuple = (;))
    unrolled_map(chain) do v
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
# Actual vapour pressure, by whichever inputs the tier carries: dewpoint (ERA5)
# → saturation vapour pressure at the dewpoint; otherwise specific humidity +
# pressure (NCEP, envelope sources).
# `hasproperty` folds to a constant for each concrete buffer tier, so the unused
# branch is eliminated.
function derive!(::Val{:actual_vapour_pressure}, buffers, ctx)
    if hasproperty(buffers, :dewpoint_temperature)
        method = ctx.vapour_pressure_method
        @inbounds for k in eachindex(buffers.actual_vapour_pressure)
            buffers.actual_vapour_pressure[k] =
                vapour_pressure(method, buffers.dewpoint_temperature[k])
        end
    else
        _vapour_pressure_from_specific_humidity!(buffers.actual_vapour_pressure,
            buffers.specific_humidity, buffers.pressure)
    end
    return nothing
end

# Used by sources that provide actual vapour pressure (e.g. WorldClim's
# `:vapr` → `actual_vapour_pressure`, or NCEP via the SH→VP derivation
# above). Converts it into the vapour pressure deficit that the relative-humidity derivations consume:
#   vapour_pressure_deficit = saturation_vapour_pressure(mean_temperature) - actual_vapour_pressure.
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
    _wind_speed!(buffers.wind_speed, buffers.eastward_wind, buffers.northward_wind)
    return nothing
end

# Power-law wind height correction: scale wind from `z_src` measurement height
# to `z_ref` reference height using a neutral-stability shear exponent α = 0.15.
# Accepts plain numbers (metres assumed) or Unitful lengths; dividing same-unit
# quantities yields a dimensionless ratio, so no stripping is needed.
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
        buffers.global_radiation, flat_terrain, buffers.days_of_year,
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

# Deep-soil temperature: the annual mean of the tier's air temperature
# (`reference_temperature` where present, else `mean_temperature`), replicated
# across that year's output days. `ctx.steps_per_year` is the output day count
# per year, so this serves monthly, daily, and hourly tiers alike.
function derive!(::Val{:deep_soil_temperature}, buffers, ctx)
    _deep_soil!(buffers.deep_soil_temperature, _air_temperature(buffers), ctx.steps_per_year)
    return nothing
end

@inline _air_temperature(b) =
    hasproperty(b, :reference_temperature) ? b.reference_temperature : b.mean_temperature

# Accumulate the native-timestep rainfall into the daily total.
function derive!(::Val{:rainfall_daily}, buffers, ctx)
    _resample!(buffers.rainfall_daily, buffers.native.rainfall,
               Accumulate(), samples_per_day(ctx.native_timestep))
    return nothing
end

# Power-law height correction of the (resampled) wind speed, in place on the
# output tier: measurement height → reference height.
function derive!(::Val{:apply_wind_shear}, buffers, ctx)
    shear = _wind_height_correction(ctx.wind_reference_height)
    @. buffers.wind_speed *= shear
    return nothing
end

# Single-value relative humidity:
#   reference_humidity = actual_vapour_pressure / saturation_vapour_pressure(reference_temperature)
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

# ---------------------------------------------------------------------------
# Solar geometry (clear-sky baseline for the radiation resample rule)
# ---------------------------------------------------------------------------

# Run solar_radiation! to populate scratch.solar.out with hourly clear-sky
# global_horizontal values for all days. Called before the shortwave
# `Disaggregate` resample, which reads the clear-sky output.
function derive!(::Val{:solar_geometry}, buffers, ctx)
    (; scratch, site) = ctx
    flat_terrain = setproperties(scratch.cloud_constants.flat_terrain_template,
        (; site.elevation,
           atmospheric_pressure = ctx.atmospheric_pressure,
           site.latitude, site.longitude))
    solar_radiation!(scratch.solar.out, scratch.solar.buffers,
        scratch.cloud_constants.solar_model;
        solar_terrain = flat_terrain,
        days  = ctx.days_of_year,
        hours = scratch.cloud_constants.hours)
    return nothing
end

# ---------------------------------------------------------------------------
# Helpers reused by derivations
# ---------------------------------------------------------------------------
#

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

# ---------------------------------------------------------------------------
# Date-range helpers
# ---------------------------------------------------------------------------

# Distinct solar-geometry days/year — sizes `scratch.solar.out`. Always the
# calendar's day count, independent of timestep.
@inline _solar_ndays_per_year(cal::Calendar) = days_per_year(cal)

@inline _minmax_env_type(::Monthly) = MonthlyMinMaxEnvironment
@inline _minmax_env_type(::Daily)   = DailyMinMaxEnvironment

@inline _days_of_year(::Monthly, nyears) = repeat(Microclimate.DEFAULT_DAYS, nyears)
@inline _days_of_year(::Daily,   nyears) = repeat(1:365, nyears)

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
# no Feb 29). The Ti slice is in units of the *native* timestep — that is how
# the loaded stack is laid out — so it is parameterised by `native_timestep`.
function _ti_range_for_dates(::Monthly, ::Timestep, years, dates_vec::Vector{Date})
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

# Daily calendar: the native stack holds `spd = samples_per_day(native_timestep)`
# steps per day (1 for MinMax, 24 for hourly, 4 for six-hourly). The returned
# day-of-year vector always has one entry per calendar day — solar geometry and
# derivations work at daily granularity regardless of sub-daily timestep.
function _ti_range_for_dates(::Daily, native::Timestep, years, dates_vec::Vector{Date})
    spd = samples_per_day(native)
    start_d, end_d = minimum(dates_vec), maximum(dates_vec)
    years_v = collect(years)
    yi_start = findfirst(==(year(start_d)), years_v)
    yi_end   = findfirst(==(year(end_d)),   years_v)
    ti_start = (yi_start - 1) * 365 * spd + (_doy_noleap(start_d) - 1) * spd + 1
    ti_end   = (yi_end   - 1) * 365 * spd + _doy_noleap(end_d) * spd
    return ti_start, ti_end, _doy_noleap.(dates_vec)
end

# One "anchor" Date per solver step — used to carry the year when building
# the output DateTime axis. Monthly: 1st of each selected month. Daily: the
# dates_vec itself (already one per calendar day).
function _step_anchor_dates(::Monthly, dates_vec::Vector{Date})
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
_step_anchor_dates(::Daily, dates_vec::Vector{Date}) = dates_vec

# Monthly forcing → independent representative days (Fortran "monthly mode");
# daily forcing → continuous run with state carrying day-to-day.
@inline _time_mode(::Monthly) = Microclimate.NonConsecutiveDayMode()
@inline _time_mode(::Daily)   = Microclimate.ConsecutiveDayMode()

