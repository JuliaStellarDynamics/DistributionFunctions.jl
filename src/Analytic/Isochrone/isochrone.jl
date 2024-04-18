#####################################
# Isochrone distribution functions (analytic)
#####################################
abstract type IsochroneDistributionFunction <: DistributionFunction end
struct IsotropicIsochrone <: IsochroneDistributionFunction
    potential::IsochronePotential # Potential model
end
struct OsipkovMerrittIsochrone <: IsochroneDistributionFunction
    ra::Float64                # Anisotropy radius
    potential::IsochronePotential # Potential model
end

"""
the Isochrone distribution function scale
"""
function dfscale(df::IsochroneDistributionFunction)
    return (df.potential.G*df.potential.M*df.potential.bc)^(-3/2)
end

include("isotropic.jl")
include("osipkovmerritt.jl")