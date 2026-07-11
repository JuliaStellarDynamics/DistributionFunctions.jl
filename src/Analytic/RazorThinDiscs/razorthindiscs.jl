
using SpecialFunctions # for gamma function
using HypergeometricFunctions

# @IMPROVE is there a better way to do this?
const IntorFloat = Union{Int64,Float64}


#####################################
# Disc distribution functions (analytic)
#####################################
abstract type DiscDF <: DiscEnergyAngularMomentumDF end
abstract type MestelPotentialDF <: DiscDF end
abstract type ZangDF <: MestelPotentialDF end
abstract type ToomrePotentialDF <: DiscDF end

const MestelPotentials = Union{MestelPotential,TaperedMestel}


# @IMPROVE, these potentials are not specific at all: but need to be either MestelPotential or TaperedMestel
struct MestelDisc{modelT<:MestelPotentials,qT<:IntorFloat} <: MestelPotentialDF
    potential::modelT # potential model
    q::qT                   # velocity dispersion parameter
    G::Float64                      # gravitational constant (not in MestelPotential or TaperedMestel, so needed here)
    isOdd::Bool
end

struct ZangDisc{modelT<:MestelPotentials,qT<:IntorFloat} <: ZangDF
    potential::modelT # potential model
    q::qT                # velocity dispersion parameter
    ν::Int64                   # inner taper
    Rin::Float64               # inner taper radius
    μ::Int64                   # outer taper
    Rout::Float64              # outer taper radius
    G::Float64                      # gravitational constant (not in MestelPotential or TaperedMestel, so needed here)  
    isOdd::Bool 
end

struct TruncatedZangDisc{modelT<:MestelPotentials,qT<:IntorFloat} <: ZangDF
    potential::modelT # potential model
    q::qT                 # velocity dispersion parameter
    ν::Int64                   # inner taper
    Rin::Float64               # inner taper radius
    μ::Int64                   # outer taper
    Rout::Float64              # outer taper radius
    Rmax::Float64              # no particles beyond Rmax
    G::Float64                      # gravitational constant (not in MestelPotential or TaperedMestel, so needed here)
    isOdd::Bool
    ξ::Float64                # Self-gravity fraction
end

struct ToomreDisc{modelT<:ToomrePotential} <: ToomrePotentialDF
    potential::modelT # potential model
    mM::Int64            # Miyamoto index
    G::Float64        # gravitational constant
    isOdd::Bool
end

struct ToomreDiscOdd{modelT<:ToomrePotential} <: ToomrePotentialDF
    potential::modelT # potential model
    mM::Int64            # Miyamoto index
    G::Float64        # gravitational constant
    isOdd::Bool
end


#####
#
#   The Mestel functions
#
#####

"""
    σMestelDistribution([R0, V0, q])

radial velocity dispersion of the tapered Mestel DF
"""
function σMestelDistribution(df::MestelPotentialDF)::Float64
    return df.potential.V0 / sqrt(df.q+1)
end

"""
    NormConstMestelDistribution([R0, V0, q])

normalization constant of the tapered Mestel DF.
"""
function NormConstMestelDistribution(df::MestelPotentialDF)::Float64
    σ = σMestelDistribution(df)
    return (df.potential.V0)^(2) / ( 2^(df.q/2+1) * (pi)^(3/2) * df.G * gamma(0.5+0.5*df.q) * (σ)^(df.q+2) * (df.potential.R0)^(df.q+1) )
end

"""
    MestelDistribution(EL::Tuple{Float64,Float64},df::MestelDisc)
Mestel distribution function.
"""
function MestelDistribution(EL::Tuple{Float64,Float64},df::MestelPotentialDF)::Float64

    E,L = EL
    σ = σMestelDistribution(df)
    C = NormConstMestelDistribution(df)

    return C * (L)^(df.q) * exp(-E / (σ^2))
end

"""
    MesteldFdE(EL::Tuple{Float64,Float64},df::MestelDisc)
Mestel DF derivative w.r.t. E.
"""
function MesteldFdE(EL::Tuple{Float64,Float64},df::MestelPotentialDF)::Float64

    σ = σMestelDistribution(df)
    return - MestelDistribution(EL,df) / (σ^2)
end

"""
    MesteldFdL(E, L[, C, q, sigma])
Mestel DF derivative w.r.t. E.
"""
function MesteldFdL(EL::Tuple{Float64,Float64},df::MestelPotentialDF)::Float64

    E,L = EL
    σ = σMestelDistribution(df)
    C = NormConstMestelDistribution(df)
    return C * df.q * (L)^(df.q-1) * exp(-E / (σ^2))
end

#####
#
#   The Miyamoto DFs (even component) for Kuzmin-Toomre disc 
#
#####

"""
    MiyamotoDistribution(EL::Tuple{Float64,Float64},df::ToomreDisc)
Miyamoto distribution function.
"""
function MiyamotoDistribution(EL::Tuple{Float64,Float64},df::ToomrePotentialDF)::Float64

    E,L = EL
    mM = df.mM

    M = df.potential.M
    E0 = energy_scale(df.potential)
    L0 = momentum_scale(df.potential)

    tE = E/E0
    tL = L/L0

    return M/(L0^2) * (
        (2mM + 3)
        * (tE)^(2mM + 2)
        * _₂F₁(-mM, -2-2mM, 1/2, tL^2/(2*tE))
        / (4 * pi^2)
    )
end

"""
    MiyamotodFdE(EL::Tuple{Float64,Float64},df::ToomreDisc)
Miyamoto DF derivative w.r.t. E.
"""
function MiyamotodFdE(EL::Tuple{Float64,Float64},df::ToomrePotentialDF)::Float64

    E,L = EL
    mM = df.mM

    M = df.potential.M
    E0 = energy_scale(df.potential)
    L0 = momentum_scale(df.potential)

    tE = E/E0
    tL = L/L0

    dtFdtE = - M/(L0^2) * (
        tE^(2mM)
        * (1 + mM)
        * (3 + 2mM)
        * ( 
            tL^2 * mM * _₂F₁(-1-2mM, 1-mM, 3/2, tL^2/(2*tE))
            - tE * _₂F₁(-mM, -2-2mM, 1/2, tL^2/(2*tE))
        )
        / (2*pi^2)
    )
    dFdE = dtFdtE / E0

    return dFdE
end

"""
    MiyamotodFdL(EL::Tuple{Float64,Float64},df::ToomreDisc)
Miyamoto DF derivative w.r.t. L.
"""
function MiyamotodFdL(EL::Tuple{Float64,Float64},df::ToomrePotentialDF)::Float64

    E,L = EL
    mM = df.mM

    M = df.potential.M
    E0 = energy_scale(df.potential)
    L0 = momentum_scale(df.potential)

    tE = E/E0
    tL = L/L0

    dtFdtL = M/(L0^2) * (
        tE^(1 + 2mM)
        * tL
        * mM
        * (1 + mM)
        * (3 + 2mM)
        * _₂F₁(-1-2mM, 1-mM, 3/2, tL^2/(2*tE))
        / pi^2
    )
    dFdL = dtFdtL / L0

    return dFdL
end




#####
#
#   The Miyamoto DFs (odd component) for Kuzmin-Toomre disc 
#
#####

"""
    MiyamotoDistributionOdd(EL::Tuple{Float64,Float64},df::ToomreDisc)
Miyamoto distribution function (odd component) .
"""
function MiyamotoDistributionOdd(EL::Tuple{Float64,Float64},df::ToomrePotentialDF)::Float64

    E,L = EL
    mM = df.mM

    M = df.potential.M
    E0 = energy_scale(df.potential)
    L0 = momentum_scale(df.potential)

    tE = E/E0
    tL = L/L0

    return M/(L0^2) * (
        (4mM + 7) * (4mM + 9)
        * (tE)^(2mM + 5/2) * tL
        * _₂F₁(-mM, -5/2-2mM, 3/2, tL^2/(2*tE))
        / (16 * pi^2)
    )
end

"""
    MiyamotodFdEOdd(EL::Tuple{Float64,Float64},df::ToomreDisc)
Miyamoto DF derivative w.r.t. E (odd component) .
"""
function MiyamotodFdEOdd(EL::Tuple{Float64,Float64},df::ToomrePotentialDF)::Float64

    E,L = EL
    mM = df.mM

    M = df.potential.M
    E0 = energy_scale(df.potential)
    L0 = momentum_scale(df.potential)

    tE = E/E0
    tL = L/L0

    dtFdtE = M/(L0^2) * (
        tE^(3/2 + 2mM)
        * tL
        * (5 + 4mM)
        * (7 + 4mM)
        * (9 + 4mM)
        * _₂F₁(-mM, -3/2-2mM, 3/2, tL^2/(2*tE))
        / (32*pi^2)
    )
    dFdE = dtFdtE / E0

    return dFdE
end

"""
    MiyamotodFdLOdd(EL::Tuple{Float64,Float64},df::ToomreDisc)
Miyamoto DF derivative w.r.t. L (odd component) .
"""
function MiyamotodFdLOdd(EL::Tuple{Float64,Float64},df::ToomrePotentialDF)::Float64

    E,L = EL
    mM = df.mM

    M = df.potential.M
    E0 = energy_scale(df.potential)
    L0 = momentum_scale(df.potential)

    tE = E/E0
    tL = L/L0

    if (mM != 0)
        dtFdtL =  M/(L0^2) * (
            tE^(2mM+3/2)
            * (7 + 4mM)
            * (9 + 4mM)
            * ( 
                tL^2 * mM * (5 + 4mM) * _₂F₁(-3/2-2mM, 1-mM, 5/2, tL^2/(2*tE))
                + 3 * tE * _₂F₁(-mM, -5/2-2mM, 3/2, tL^2/(2*tE))
            )
            / (48*pi^2)
        )
    else
        dtFdtL =  63 * M/(L0^2) * (
            tE^(2mM+5/2)
            / (16*pi^2)
        )
    end
    dFdL = dtFdtL / L0

    return dFdL
end


include("mestel.jl")
include("zang.jl")
include("truncatedzang.jl")
include("miyamoto.jl")