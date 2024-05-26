#####################################
# Plummer distribution functions (analytic)
#####################################
struct IsotropicPlummer{modelT<:PlummerPotential} <: EnergyOnlyDistributionFunction
    potential::modelT # Potential model
end
struct OsipkovMerrittPlummerEL{modelT<:PlummerPotential} <: SphericalEnergyAngularMomentumDistributionFunction
    ra::Float64                # Anisotropy radius
    potential::modelT # Potential model
end
struct OsipkovMerrittPlummerJL{modelT<:PlummerPotential} <: SphericalActionDistributionFunction
    ra::Float64                # Anisotropy radius
    potential::modelT # Potential model
end

# create a type union for all Plummer distribution functions
const PlummerDistributionFunction = Union{IsotropicPlummer,OsipkovMerrittPlummerEL,OsipkovMerrittPlummerJL}

# unify the osipkov-merritt methods
const OsipkovMerrittPlummer = Union{OsipkovMerrittPlummerEL,OsipkovMerrittPlummerJL}

"""
the Plummer distribution function scale
"""
function dfscale(df::PlummerDistributionFunction)
    return (df.potential.G*df.potential.M*df.potential.bc)^(-3/2)
end

include("isotropic.jl")
include("osipkovmerrit.jl")