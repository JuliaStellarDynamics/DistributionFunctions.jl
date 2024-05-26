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

# unify all energy angular momentum types
const EnergyAngularMomentumDistributionFunction = Union{EnergyOnlyDistributionFunction,SphericalEnergyAngularMomentumDistributionFunction,DiscEnergyAngularMomentumDistributionFunction}

# unify all action types
const ActionDistributionFunction = Union{SphericalActionDistributionFunction,DiscActionDistributionFunction}


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
    gradient(EL::Tuple{Float64,Float64}, distributionfunction::EnergyAngularMomentumDistributionFunction)

Angular momentum derivative of a given distribution function `distributionfunction` for a given `E`,`L`.
"""
function gradient(EL::Tuple{Float64,Float64}, df::EnergyAngularMomentumDistributionFunction)

    DFDEval = DFDE(EL,df)
    DFDLval = DFDL(EL,df)

    return DFDEval,DFDLval
    
end

"""
    gradient(JL::Tuple{Float64,Float64}, distributionfunction::ActionDistributionFunction)

Angular momentum derivative of a given distribution function `distributionfunction` for a given `Jr`,`L`.

"""
function gradient(EL::Tuple{Float64,Float64}, df::ActionDistributionFunction; Ω1::Float64=-1)

    DFDEval = DFDE(EL,df)
    DFDLval = DFDL(EL,df)

    if (Ω1 == -1)
        a,e = ae_from_EL(E,L,df.potential)
        Ω1,_ = frequencies_from_ae(a,e,df.potential)
    end

    # now convert using the Jacobian dE/dJ = Ω1
    return Ω1 * DFDEval,DFDLval
    
end


#####################################
# Include analytic distribution functions
#####################################
include("../Analytic/Spheres/Isochrone/isochrone.jl")
include("../Analytic/Spheres/Plummer/plummer.jl")
include("../Analytic/RazorThinDiscs/razorthindiscs.jl")