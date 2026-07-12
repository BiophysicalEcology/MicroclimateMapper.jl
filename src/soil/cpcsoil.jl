# CPCSoil — CPC Global Monthly Soil Moisture.
#
# A single NetCDF with the `:soilw` layer holding 12 monthly steps on a 0–360°
# longitude grid. Declared like any source providing one canonical variable
# (soil moisture): the SingleFileBands loader reads the file, the Longitude360
# convention handles the 0–360° grid, and the 12 months tile across the run's
# years. Loaded as a prescribed forcing via `_load_prescribed`.

loader(::Type{CPCSoil}) = SingleFileBands()
longitude_convention(::Type{CPCSoil}) = Longitude360()

_extra_getraster_kwargs(::Type{CPCSoil}) = (; period = "1991-2020")

# Native `:soilw` is mm of soil water in a 1 m column; mm/m converts that storage
# to the dimensionless volumetric fraction m³/m³ (1 mm / 1 m = 10⁻³).
variables(::Type{CPCSoil}) = (Variable(SoilMoisture(), :soilw, u"mm/m"),)
