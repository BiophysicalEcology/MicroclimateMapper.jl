# Soil texture fetch, shared between SoilGrids and SLGA (same 4 layers, same
# 6 depth bins). `TextureVariable` mirrors `weather.jl`'s `WeatherVariable`,
# kept separate to avoid coupling to the weather-derivation machinery.

struct TextureVariable{Name, Field, U, T}
    unit::U
    transform::T
end
TextureVariable(name::Symbol, field::Symbol, unit = 1, transform = identity) =
    TextureVariable{name, field, typeof(unit), typeof(transform)}(unit, transform)

canonical_name(::TextureVariable{Name}) where {Name} = Name
native_field(::TextureVariable{<:Any, Field}) where {Field} = Field

texture_variables(::Type{SoilGrids}) = (
    TextureVariable(:bulk_density, :bdod, u"Mg/m^3", raw -> raw / 100),  # cg/cm^3 -> Mg/m^3
    TextureVariable(:clay, :clay),  # g/100g already == %
    TextureVariable(:silt, :silt),
    TextureVariable(:sand, :sand),
)

# Reduce each depth-bin raster to one value (single uniform profile, not per-pixel).
function _texture_values_from_paths(paths, area::Extent, var::TextureVariable)
    map(paths) do path
        window = read(crop(Raster(path; lazy = true); to = area, touches = true))
        var.transform(mean(skipmissing(window))) * var.unit
    end
end

function _load_soil_texture_native(::Type{SoilGrids}, area::Extent; quantile = "mean")
    vars = texture_variables(SoilGrids)
    depth_bins = collect(depths(SoilGrids))  # getraster's depth::AbstractArray dispatch needs a Vector, not a Tuple
    values = map(vars) do var
        paths = getraster(SoilGrids, native_field(var); depth = depth_bins, quantile)
        _texture_values_from_paths(paths, area, var)
    end
    return NamedTuple{map(canonical_name, vars)}(values)
end

# ISRIC's point REST API -- a different endpoint from the raster VRTs.
function _load_soil_texture_native(::Type{SoilGrids}, lon::Real, lat::Real; quantile = "mean")
    vars = texture_variables(SoilGrids)
    depth_bins = depths(SoilGrids)
    values = map(vars) do var
        map(depth_bins) do depth
            point = PointDataSources.getpoint(SoilGrids, native_field(var); lon, lat, depth, quantile)
            _check_soilgrids_units(var, point.units)
            var.transform(point.value) * var.unit
        end
    end
    return NamedTuple{map(canonical_name, vars)}(values)
end

function _check_soilgrids_units(::TextureVariable{Name}, units::AbstractString) where {Name}
    Name === :bulk_density && !occursin("cg", units) && !occursin("g/cm", units) &&
        @warn "SoilGrids point API returned unexpected units \"$units\" for bulk density; expected cg/cm^3."
    return nothing
end
