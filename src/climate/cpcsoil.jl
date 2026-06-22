# cpcsoil.jl
#
# Companion loader for CPC Global Monthly Soil Moisture (CPCSoil).
# Intended for use alongside CRUCL2. Pass the result as `data.soil_moisture`
# on a MicroRasterProblem or MicroVectorProblem.

# 1 mm water column per 1 m soil depth = 10⁻³ m³/m³ (native CPCSoil unit).
const _CPCSOIL_MM_TO_VOLUMETRIC = 1e-3

weather_variables(::Type{CPCSoil}) = (
    WeatherVariable(:soil_moisture, :soilw, 1, raw -> raw * _CPCSOIL_MM_TO_VOLUMETRIC),
)

"""
    load_cpcsoil(area, years; period="1991-2020") -> Raster

Load CPC Global Monthly Soil Moisture over `area`, tiled to `12 * length(years)`
time steps matching a CRUCL2 weather stack. Returns a `(X, Y, Ti)` Raster in
m³/m³. `period` is `"1991-2020"` (default) or `"1981-2010"`.
"""
_load_soil_moisture(::Type{CPCSoil}, area, years) = load_cpcsoil(area, years)

function load_cpcsoil(area::Extent, years; period = "1991-2020")
    path   = getraster(CPCSoil; period)
    nyears = length(years)

    # CPCSoil uses 0–360° longitude. Shift the crop area into 0-360,
    # then shift the loaded coordinates back to -180-180.
    needs_shift = area.X[1] < 0 || area.X[2] < 0
    load_area   = needs_shift ?
        Extent(X = (mod(area.X[1], 360.0), mod(area.X[2], 360.0)), Y = area.Y) :
        area

    raw = Rasters.replace_missing(
        read(crop(Raster(path; name = :soilw, lazy = true); to = load_area, touches = true)),
        NaN)

    if needs_shift
        x_lk = lookup(raw, X)
        new_x = Sampled(collect(x_lk) .- 360.0;
            order = order(x_lk), span = span(x_lk), sampling = sampling(x_lk))
        raw = Raster(parent(raw), (X(new_x), dims(raw, Y), dims(raw, Ti)); crs = crs(raw))
    end

    # Apply the declared transform (native mm/m → m³/m³), matching how
    # _read_one_variable! would apply WeatherVariable.transform for a native source.
    sm_transform = only(weather_variables(CPCSoil)).transform
    data_m3 = sm_transform.(parent(raw))
    tiled   = nyears == 1 ? data_m3 : repeat(data_m3; outer = (1, 1, nyears))

    return Raster(tiled, (dims(raw)[1:2]..., Ti(1:(12 * nyears))); crs = crs(raw))
end
