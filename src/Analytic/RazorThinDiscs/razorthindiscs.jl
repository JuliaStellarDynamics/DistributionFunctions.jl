
using SpecialFunctions # for gamma function

# @IMPROVE is there a better way to do this?
IntorFloat = Union{Int64,Float64}


#####################################
# Disc distribution functions (analytic)
#####################################
abstract type DiscDistributionFunction <: DistributionFunction end
abstract type MestelPotentialDistributionFunction <: DiscDistributionFunction end
abstract type ZangDistributionFunction <: MestelPotentialDistributionFunction end

# @IMPROVE, these potentials are not specific at all: but need to be either MestelPotential or TaperedMestel
struct MestelDisc <: MestelPotentialDistributionFunction
    potential::CentralCorePotential # potential model
    q::IntorFloat                   # dispersion parameter
    ξDF::Float64                    # self-gravity fraction
    G::Float64                      # gravitational constant (not in MestelPotential or TaperedMestel, so needed here)
end

struct ZangDisc <: ZangDistributionFunction
    potential::CentralCorePotential # potential model
    q::IntorFloat                # dispersion parameter
    ν::Int64                   # inner taper
    Rin::Float64               # inner taper radius
    μ::Int64                   # outer taper
    Rout::Float64              # outer taper radius
    ξDF::Float64               # self-gravity fraction 
    G::Float64                      # gravitational constant (not in MestelPotential or TaperedMestel, so needed here)   
end

struct TruncatedZangDisc <: ZangDistributionFunction
    potential::CentralCorePotential # potential model
    q::IntorFloat                 # dispersion parameter
    ν::Int64                   # inner taper
    Rin::Float64               # inner taper radius
    μ::Int64                   # outer taper
    Rout::Float64              # outer taper radius
    Rmax::Float64              # no particles beyond Rmax
    ξDF::Float64               # self-gravity fraction  
    G::Float64                      # gravitational constant (not in MestelPotential or TaperedMestel, so needed here)
end

"""
    σMestelDF([R0, V0, q])

radial velocity dispersion of the tapered Mestel DF
"""
function σMestelDF(df::MestelPotentialDistributionFunction)::Float64
    return df.potential.V0 / sqrt(df.q+1)
end

"""
    NormConstMestelDF([R0, V0, q])

normalization constant of the tapered Mestel DF.
"""
function NormConstMestelDF(df::MestelPotentialDistributionFunction)::Float64
    σ = σMestelDF(df)
    return (df.potential.V0)^(2) / ( 2^(df.q/2+1) * (pi)^(3/2) * df.G * gamma(0.5+0.5*df.q) * (σ)^(df.q+2) * (df.potential.R0)^(df.q+1) )
end

"""
    MestelDF(EL::Tuple{Float64,Float64},df::MestelDisc)
Mestel distribution function.
"""
function MestelDF(EL::Tuple{Float64,Float64},df::MestelPotentialDistributionFunction)::Float64

    E,L = EL
    σ = σMestelDF(df)
    C = NormConstMestelDF(df)

    return C * (L)^(df.q) * exp(-E / (σ^2))
end

"""
    MesteldFdE(EL::Tuple{Float64,Float64},df::MestelDisc)
Mestel DF derivative w.r.t. E.
"""
function MesteldFdE(EL::Tuple{Float64,Float64},df::MestelPotentialDistributionFunction)::Float64

    σ = σMestelDF(df)
    return - F(EL,df) / (σ^2)
end

"""
    MesteldFdL(E, L[, C, q, sigma])
Mestel DF derivative w.r.t. E.
"""
function MesteldFdL(EL::Tuple{Float64,Float64},df::MestelPotentialDistributionFunction)::Float64

    E,L = EL
    σ = σMestelDF(df)
    C = NormConstMestelDF(df)
    return C * df.q * (L)^(df.q-1) * exp(-E / (σ^2))
end

include("mestel.jl")
include("zang.jl")
include("truncatedzang.jl")
#include("miyamoto.jl")