#####################################
# Abstract types
#####################################
"""
Abstract type for distribution functions
"""
abstract type DistributionFunction end

"""
Abstract type for energy-only distribution functions
"""
abstract type EnergyOnlyDistributionFunction <: DistributionFunction end


#####################################
# Generic functions (not methods)
#####################################
"""
    F(EL::Tuple{Float64,Float64}, distributionfunction::DistributionFunction)

Distribution function `distributionfunction` for a given `E`,`L`.
"""
function F(EL::Tuple{Float64,Float64}, df::DistributionFunction)
    # ... [model specific implementation] ...
end

"""
    ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},res::Resonance, distributionfunction::DistributionFunction)

Distribution function derivative for `distributionfunction` for a given `E`,`L`.
"""
function ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::DistributionFunction)
        # ... [model specific implementation] ...
end


#####################################
# Include analytic distribution functions
#####################################
include("../Analytic/Isochrone/isochrone.jl")
include("../Analytic/Plummer/plummer.jl")
include("../Analytic/RazorThinDiscs/razorthindiscs.jl")