module MicroclimateMapperMakieExt

using Makie
using MicroclimateMapper
using MicroclimateMapper: MicroMapCache
using Rasters
using Rasters: AbstractRaster, X, Y, Ti, hasdim, dims, rebuild
using Unitful: Quantity, unit, ustrip

# Spatial preview of a `MicroMapCache`'s **inputs**: the data that's been
# loaded and pre-resolved by `init(problem)` but not yet consumed by
# `solve!`. Each layer is rendered via `Rasters.rplot`. A shared time
# slider drives the time-varying layers (weather, canonical overrides);
# static layers (terrain, surface properties) ignore it.
#
# The plot does no simulation work — it only reads `cache` and reacts to
# slider events.
#
# Pass `layers` to restrict what's shown to a subset of the input layer
# names (see `MicroclimateMapper.input_layer_names(cache)`). Extra keyword
# arguments are forwarded to `Rasters.rplot`.
function Makie.plot(cache::MicroMapCache; layers = nothing, rplot_kw...)
    fig = Figure()
    _plot_inputs!(fig, cache, layers; rplot_kw...)
    return fig
end

function _plot_inputs!(fig, cache, requested_layers; rplot_kw...)
    all_layers = _collect_input_layers(cache)
    selected   = _select_layers(all_layers, requested_layers)
    isempty(selected) && return fig

    # Build the time slider only if any selected layer varies in time.
    nt = _max_time_length(selected)
    plots_layout = fig[1, 1] = GridLayout()
    sl_time = if nt > 0
        slider_layout = fig[2, 1] = GridLayout(tellwidth = false)
        sl = Slider(slider_layout[1, 2], range = 1:nt, startvalue = 1)
        Label(slider_layout[1, 1], "time")
        Label(slider_layout[1, 3], lift(t -> string(" step ", t), sl.value))
        sl
    else
        nothing
    end

    nrows, ncols = _balance_grid(length(selected))
    for (i, (name, raster)) in enumerate(pairs(selected))
        slice_obs = _input_slice_observable(raster, sl_time)
        u = _layer_unit(raster)
        # Makie's UnitfulConversion can't handle affine units (e.g. °C),
        # so strip units and fold them into the title.
        plot_obs = u === nothing ? slice_obs :
                   lift(r -> rebuild(r, ustrip.(u, parent(r))), slice_obs)
        title = u === nothing ? string(name) : "$name ($u)"
        row, col = fldmod1(i, ncols)
        Rasters.rplot(plots_layout[row, col], plot_obs; title, rplot_kw...)
    end

    return fig
end

# Gather every input layer the cache holds — terrain (skipping the
# `horizon_angles` SVector field which isn't a heatmap), pre-resolved
# surface grids, the native weather stack, and any user `data` overrides.
function _collect_input_layers(cache)
    pairs = Pair{Symbol, Any}[]

    # Terrain — only the per-pixel scalar fields. `horizon_angles` is an
    # SVector per pixel and `latitude`/`longitude` are just X/Y broadcast;
    # neither is useful as a 2-D map here.
    terrain = cache.terrain
    if terrain !== nothing
        for name in (:elevation, :slope, :aspect, :atmospheric_pressure)
            hasproperty(terrain, name) || continue
            push!(pairs, name => getproperty(terrain, name))
        end
    end

    cache.albedo_grid    === nothing || push!(pairs, :surface_albedo   => cache.albedo_grid)
    cache.roughness_grid === nothing || push!(pairs, :roughness_height => cache.roughness_grid)

    # Native weather stack — use canonical variable names as labels (e.g.
    # `:maximum_temperature` rather than the source's `:tmax`). The mapping
    # comes from `weather_variables(weather_source)`. Dedup by native field
    # so a raster that's mapped to multiple canonical names appears once.
    weather = cache.weather
    if weather !== nothing
        source = cache.problem.model.weather_source
        if source === nothing
            for name in propertynames(weather)
                push!(pairs, name => getproperty(weather, name))
            end
        else
            seen = Set{Symbol}()
            for var in MicroclimateMapper.weather_variables(source)
                field = MicroclimateMapper.native_field(var)
                field in seen && continue
                push!(seen, field)
                canonical = MicroclimateMapper.canonical_name(var)
                push!(pairs, canonical => getproperty(weather, field))
            end
        end
    end

    # User canonical-variable overrides from `problem.data` (already
    # resampled to the run template at `init` time).
    for name in keys(cache.canonical_overrides)
        push!(pairs, name => cache.canonical_overrides[name])
    end

    return NamedTuple(pairs)
end

_select_layers(all_layers, ::Nothing) = all_layers
_select_layers(all_layers, requested) =
    NamedTuple{Tuple(requested)}(map(k -> all_layers[k], requested))

# Wrap an input raster as a time-sliceable Observable. Static (no `Ti`)
# rasters are returned wrapped in a constant Observable so the same
# downstream pipeline handles both.
function _input_slice_observable(raster, sl_time)
    if sl_time !== nothing && hasdim(raster, Ti)
        return lift(t -> view(raster, Ti(t)), sl_time.value)
    else
        return Observable(raster)
    end
end

_max_time_length(layers) =
    reduce(max, (length(dims(r, Ti)) for r in values(layers) if hasdim(r, Ti)); init = 0)

# Native unit (or `nothing`) of a Unitful-valued raster.
_layer_unit(r::AbstractRaster{<:Quantity}) = unit(eltype(r))
_layer_unit(::AbstractRaster) = nothing

# Rough square layout (favouring more columns than rows).
function _balance_grid(n::Int)
    n == 0 && return (1, 1)
    ncols = ceil(Int, sqrt(n))
    nrows = ceil(Int, n / ncols)
    return nrows, ncols
end

end
