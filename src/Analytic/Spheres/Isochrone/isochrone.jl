#####################################
# Isochrone distribution functions (analytic)
#####################################

struct IsotropicIsochrone{modelT<:IsochronePotential} <: ErgodicDF
    potential::modelT # Potential model
end
struct OsipkovMerrittIsochroneEL{modelT<:IsochronePotential} <: SphericalEnergyAngularMomentumDF
    ra::Float64                # Anisotropy radius
    potential::modelT # Potential model
end
struct OsipkovMerrittIsochroneJL{modelT<:IsochronePotential} <: SphericalEnergyAngularMomentumDF
    ra::Float64                # Anisotropy radius
    potential::modelT # Potential model
end

# create a type union for all isochrone distribution functions
const IsochroneDF = Union{IsotropicIsochrone,OsipkovMerrittIsochroneEL,OsipkovMerrittIsochroneJL}

# unify the osipkov-merritt models
const OsipkovMerrittIsochrone = Union{OsipkovMerrittIsochroneEL,OsipkovMerrittIsochroneJL}

"""
the Isochrone distribution function scale
"""
function dfscale(df::IsochroneDF)
    return (df.potential.G*df.potential.M*df.potential.bc)^(-3/2)
end

include("isotropic.jl")
include("osipkovmerritt.jl")