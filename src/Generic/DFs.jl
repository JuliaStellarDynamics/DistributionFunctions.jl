#####################################
# Abstract types
#####################################
"""
Abstract type for distribution functions
"""
abstract type DistributionFunction end

"""
Geometries
"""
abstract type SphericalDistributionFunction <: DistributionFunction end
abstract type RazorThinDiscDistributionFunction <: DistributionFunction end

"""
Defining Coordinates
"""
abstract type EnergyOnlyDistributionFunction <: SphericalDistributionFunction end
abstract type SphericalEnergyAngularMomentumDistributionFunction <: SphericalDistributionFunction end
abstract type SphericalActionDistributionFunction <: SphericalDistributionFunction end

abstract type DiscEnergyAngularMomentumDistributionFunction <: RazorThinDiscDistributionFunction end
abstract type DiscActionDistributionFunction <: RazorThinDiscDistributionFunction end

#####################################
# Generic functions
#####################################
"""
    F(EL::Tuple{Float64,Float64}, distributionfunction::DistributionFunction)

Distribution function `distributionfunction` for a given `E`,`L`.
"""
function Distribution(EL::Tuple{Float64,Float64}, df::DistributionFunction)
    # ... [model specific implementation] ...
end

"""
    DFDE(EL::Tuple{Float64,Float64}, distributionfunction::DistributionFunction)

Energy derivative of a given distribution function `distributionfunction` for a given `E`,`L`.
"""
function DFDE(EL::Tuple{Float64,Float64}, df::DistributionFunction)
    # ... [model specific implementation] ...
end

"""
    DFDL(EL::Tuple{Float64,Float64}, distributionfunction::DistributionFunction)

Angular momentum derivative of a given distribution function `distributionfunction` for a given `E`,`L`.
"""
function DFDL(EL::Tuple{Float64,Float64}, df::DistributionFunction)
    # ... [model specific implementation] ...
end

"""
    gradient(EL::Tuple{Float64,Float64}, distributionfunction::DistributionFunction)

Angular momentum derivative of a given distribution function `distributionfunction` for a given `E`,`L`.
"""
function gradient(EL::Tuple{Float64,Float64}, df::DistributionFunction)

    DFDEval = DFDE(EL,df)
    DFDLval = DFDL(EL,df)

    return DFDEval,DFDLval
    
end

"""
    ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},res::Resonance, distributionfunction::DistributionFunction)

Distribution function derivative for `distributionfunction` for a given `E`,`L`.
"""
function ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::DistributionFunction)

        DFDEval = DFDE(EL,df)
        DFDLval = DFDL(EL,df)
        Ω1,Ω2   = ΩΩ
        n1,n2   = resonance.number[1],resonance.number[2]
        ndotΩ   = n1*Ω1 + n2*Ω2

        return DFDEval * ndotΩ + DFDLval * n2

    end


#####################################
# Include analytic distribution functions
#####################################
include("../Analytic/Isochrone/isochrone.jl")
include("../Analytic/Plummer/plummer.jl")
include("../Analytic/RazorThinDiscs/razorthindiscs.jl")