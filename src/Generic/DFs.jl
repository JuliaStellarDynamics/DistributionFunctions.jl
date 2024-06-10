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
abstract type SphericalDF <: DistributionFunction end
abstract type RazorThinDiscDF <: DistributionFunction end

"""
Defining Coordinates
"""
abstract type ErgodicDF <: SphericalDF end
abstract type SphericalEnergyAngularMomentumDF <: SphericalDF end
abstract type SphericalActionDF <: SphericalDF end

abstract type DiscEnergyAngularMomentumDF <: RazorThinDiscDF end
abstract type DiscActionDF <: RazorThinDiscDF end

# unify all energy angular momentum types
const EnergyAngularMomentumDF = Union{ErgodicDF,SphericalEnergyAngularMomentumDF,DiscEnergyAngularMomentumDF}

# shorted the calls:
#SphereELDF

# unify all action types
const ActionDF = Union{SphericalActionDF,DiscActionDF}


#####################################
# Generic functions
#####################################
"""
    F(EL::Tuple{Float64,Float64}, distributionfunction::DistributionFunction)

Distribution function `distributionfunction` for a given `E`,`L`.
"""
function DistributionFunction(EL::Tuple{Float64,Float64}, df::DistributionFunction)
    # ... [model specific implementation] ...
end


"""
    gradient(EL::Tuple{Float64,Float64}, distributionfunction::EnergyAngularMomentumDistributionFunction)

Angular momentum derivative of a given distribution function `distributionfunction` for a given `E`,`L`.
"""
function gradient(EL::Tuple{Float64,Float64}, df::EnergyAngularMomentumDF)

    DFDEval = DFDE(EL,df)
    DFDLval = DFDL(EL,df)

    return DFDEval,DFDLval
    
end

"""
    gradient(JL::Tuple{Float64,Float64}, distributionfunction::ActionDistributionFunction)

Angular momentum derivative of a given distribution function `distributionfunction` for a given `Jr`,`L`.

"""
function gradient(JL::Tuple{Float64,Float64}, df::ActionDF)

    DFDJval = DFDJ(JL,df)
    DFDLval = DFDL(JL,df)

    return DFDJval,DFDLval
    
end


#####################################
# Include analytic distribution functions
#####################################
include("../Analytic/Spheres/Isochrone/isochrone.jl")
include("../Analytic/Spheres/Plummer/plummer.jl")
include("../Analytic/RazorThinDiscs/razorthindiscs.jl")