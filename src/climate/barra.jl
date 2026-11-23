# BARRA bindings: the Bureau of Meteorology's Australian regional reanalysis
# (BARRA-R2 ~11 km, BARRA-C2 ~4 km). One NetCDF per variable per month.
#
# `BARRA{P,D}` and `BARRA{P,D,Hour}` are hourly. `BARRA{P,D,Day}` is daily
# min/max. Static layers such as `:orog` come from `BARRA{P,D,Static}`.

const _BARRADay = BARRA{<:RasterDataSources.BARRAProduct,<:RasterDataSources.BARRADomain,Day}

weather_calendar(::Type{<:BARRA}) = Daily()
native_timestep(::Type{<:BARRA}) = Hourly()
native_timestep(::Type{<:_BARRADay}) = MinMax()
loader(::Type{<:BARRA}) = MonthlyTimeSeries()
static_source(::Type{<:BARRA{P,D}}) where {P<:RasterDataSources.BARRAProduct,D<:RasterDataSources.BARRADomain} =
    BARRA{P,D,RasterDataSources.Static}

function variables(::Type{<:BARRA})
    (
        Variable(Reference(Temperature()),      :tas,     u"K"),
        Variable(WindSpeed(),                   :sfcWind, u"m/s"),
        Variable(Reference(RelativeHumidity()), :hurs,    1, percent_to_fraction),
        Variable(SeaLevelPressure(),            :psl,     u"Pa"),
        Variable(GlobalRadiation(),             :rsds,    u"W/m^2"),
        Variable(LongwaveRadiation(),           :rlds,    u"W/m^2"),
        # `:pr` (kg/m²/s) × 3600 s → kg/m² per hour.
        Variable(Rainfall(),                    :pr,      u"kg/m^2", raw -> raw * 3600.0),
        Variable(Elevation(),                   :orog,    u"m"),
    )
end

function variables(::Type{<:_BARRADay})
    (
        Variable(Temperature(Maximum()), :tasmax,  u"K"),
        Variable(Temperature(Minimum()), :tasmin,  u"K"),
        # Daily mean relative humidity; the min/max follow from vapour pressure.
        Variable(RelativeHumidity(),     :hurs,    1, percent_to_fraction),
        Variable(WindSpeed(),            :sfcWind, u"m/s"),
        Variable(GlobalRadiation(),      :rsds,    u"W/m^2"),
        # `:pr` (kg/m²/s) × 86400 s → kg/m² per day.
        Variable(Rainfall(),             :pr,      u"kg/m^2", raw -> raw * 86400.0),
        Variable(Elevation(),            :orog,    u"m"),
    )
end

# BARRA's own `:orog` grid as the DEM — same role as CRUCL2's `:elv`. 0.3°
# buffer matches SRTM's, comfortably wider than a BARRA-R2 cell (~11 km).
function _load_dem(T::Type{<:BARRA}, area::Extent)
    buffered = Extents.buffer(area, (X = 0.3, Y = 0.3))
    static = static_source(T)
    orog = _read_layer(static, :orog, getraster(static, :orog), buffered)
    return Rasters.replace_missing(orog, 0)
end
