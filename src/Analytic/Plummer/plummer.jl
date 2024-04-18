#####################################
# Plummer distribution functions (analytic)
#####################################
abstract type PlummerDistributionFunction <: DistributionFunction end
struct IsotropicPlummer <: PlummerDistributionFunction
    potential::PlummerPotential # Potential model
end
struct OsipkovMerrittPlummer <: PlummerDistributionFunction
    ra::Float64                # Anisotropy radius
    potential::PlummerPotential # Potential model
end

"""
the Plummer distribution function scale
"""
function dfscale(df::PlummerDistributionFunction)
    return (df.potential.G*df.potential.M*df.potential.bc)^(-3/2)
end

include("isotropic.jl")
include("osipkovmerrit.jl")