#####################################
# Kuzmin-Kuzutov distribution functions (analytic)
#####################################
struct OblateStackel{modelT<:KuzminKuzutovPotential} <: OblateEnergyAngularMomentumDF
    potential::modelT # Potential model
end

"""
the Kuzmin-Kuzutov distribution function scale
"""
function dfscale(df::PlummerDF)
    # tc = c/(a+c)
    # ta = a/(a+c)
    return df.potential.M/(df.potential.G*df.potential.M*(df.potential.a+df.potential.c))^(3/2)*(df.potential.tc^2)/(2^(3/2)*pi^3*df.potential.ta)
end
