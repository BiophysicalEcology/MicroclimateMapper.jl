# Ported from NicheMapR's pedotransfer.R (Cosby et al. 1984 Tables 4/5;
# Campbell 1985, Soil Physics with BASIC).
#
# R's air-entry potential is negative; campbell.jl:180 negates a positive
# magnitude at point of use. Each model returns a signed value privately;
# `pedotransfer` applies `abs` once, in one place.

abstract type PedotransferModel end
struct CosbyUnivariate <: PedotransferModel end
struct CosbyMultivariate <: PedotransferModel end
struct Campbell1985 <: PedotransferModel end

@inline function _pedotransfer_raw(::CosbyUnivariate, clay_percent, silt_percent, sand_percent, bulk_density)
    campbell_b_intermediate = clay_percent * 0.1590 + 2.910
    campbell_b = log10(10^campbell_b_intermediate / 10.2)
    air_entry_intermediate = sand_percent * -0.0131 + 1.880
    air_entry_potential_signed = 10^air_entry_intermediate / 10.2 * -1
    conductivity_intermediate = sand_percent * 0.0153 - 0.884
    saturated_conductivity = 10^conductivity_intermediate * 0.0007196666
    return (; campbell_b, air_entry_potential_signed, saturated_conductivity)
end

@inline function _pedotransfer_raw(::CosbyMultivariate, clay_percent, silt_percent, sand_percent, bulk_density)
    campbell_b_intermediate = clay_percent * 0.1570 + sand_percent * (-0.0030) + 3.10
    campbell_b = log10(10^campbell_b_intermediate / 10.2)
    air_entry_intermediate = sand_percent * (-0.0095) + silt_percent * 0.0063 + 1.54
    air_entry_potential_signed = 10^air_entry_intermediate / 10.2 * -1
    conductivity_intermediate = sand_percent * 0.0126 + clay_percent * (-0.0064) - 0.60
    saturated_conductivity = 10^conductivity_intermediate * 0.0007196666
    return (; campbell_b, air_entry_potential_signed, saturated_conductivity)
end

@inline function _pedotransfer_raw(::Campbell1985, clay_percent, silt_percent, sand_percent, bulk_density)
    clay_diameter, silt_diameter, sand_diameter = 0.001, 0.026, 1.05  # mm
    clay_fraction, silt_fraction, sand_fraction = clay_percent / 100, silt_percent / 100, sand_percent / 100
    log_mean = clay_fraction * log(clay_diameter) + sand_fraction * log(sand_diameter) + silt_fraction * log(silt_diameter)
    log_variance = clay_fraction * log(clay_diameter)^2 + sand_fraction * log(sand_diameter)^2 +
                   silt_fraction * log(silt_diameter)^2 - log_mean^2
    geometric_mean_diameter = exp(log_mean)
    geometric_std_diameter = exp(sqrt(log_variance))
    reference_air_entry_potential = -0.5 * geometric_mean_diameter^(-1 / 2)  # at bulk_density = 1.3 Mg/m^3
    campbell_b = -2 * reference_air_entry_potential + 0.2 * geometric_std_diameter
    bulk_density_value = ustrip(u"Mg/m^3", bulk_density)
    air_entry_potential_signed = reference_air_entry_potential * (bulk_density_value / 1.3)^(0.67 * campbell_b)
    saturated_conductivity = 0.004 * (1.3 / bulk_density_value)^(1.3 * campbell_b) *
                             exp(-6.9 * clay_fraction - 3.7 * silt_fraction)
    return (; campbell_b, air_entry_potential_signed, saturated_conductivity)
end

"""
    pedotransfer(model::PedotransferModel, clay_percent, silt_percent, sand_percent, bulk_density)
        -> (; campbell_b, air_entry_potential, saturated_conductivity)

Soil texture to Campbell hydraulic parameters via `model`: `CosbyUnivariate`
(Cosby et al. 1984, Table 5), `CosbyMultivariate` (Table 4), or
`Campbell1985` (particle-size theory; the only one using `bulk_density`).
`air_entry_potential` (`u"J/kg"`, positive) and `saturated_conductivity`
(`u"kg*s/m^3"`).
"""
@inline function pedotransfer(model::PedotransferModel, clay_percent, silt_percent, sand_percent, bulk_density)
    raw = _pedotransfer_raw(model, clay_percent, silt_percent, sand_percent, bulk_density)
    return (;
        campbell_b = raw.campbell_b,
        air_entry_potential = abs(raw.air_entry_potential_signed) * u"J/kg",
        saturated_conductivity = raw.saturated_conductivity * u"kg*s/m^3",
    )
end

"""
    field_capacity_and_wilting_point(clay_percent, silt_percent) -> (; field_capacity, wilting_point)

Rab et al. (2011) field capacity / wilting point (`u"m^3/m^3"`), model-
independent. Not consumed by `SoilProfile`; returned for reference.
"""
@inline function field_capacity_and_wilting_point(clay_percent, silt_percent)
    field_capacity = (7.561 + 1.176 * clay_percent - 0.009843 * clay_percent^2 + 0.2132 * silt_percent) / 100
    wilting_point = (-1.304 + 1.117 * clay_percent - 0.009309 * clay_percent^2) / 100
    return (; field_capacity = field_capacity * u"m^3/m^3", wilting_point = wilting_point * u"m^3/m^3")
end
