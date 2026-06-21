# cpcsoil.jl
#
# Companion loader for CPC Global Monthly Soil Moisture (CPCSoil).
# Intended for use alongside CRUCL2. Pass the result as `data.soil_moisture`
# on a MicroRasterProblem or MicroVectorProblem.

"""
    load_cpcsoil(area, years; period="1991-2020") -> Raster

Load CPC Global Monthly Soil Moisture over `area`, tiled to `12 * length(years)`
time steps matching a CRUCL2 weather stack. Returns a `(X, Y, Ti)` Raster in
m³/m³ (converted from native mm storage). `period` is `"1991-2020"` (default)
or `"1981-2010"`.
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

    # Native unit: mm of soil water per 1 m column. Divide by 1 m to get
    # dimensionless volumetric fraction (m³/m³): 1 mm / 1 m = 10⁻³ m³/m³.
    data_m3 = ustrip.(u"m^3/m^3", parent(raw) .* u"mm" ./ 1u"m")
    tiled   = nyears == 1 ? data_m3 : repeat(data_m3; outer = (1, 1, nyears))

    return Raster(tiled, (dims(raw)[1:2]..., Ti(1:(12 * nyears))); crs = crs(raw))
end
