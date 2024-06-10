

"""
MestelDisc([potential])

Mestel disc distribution function.
"""
function MestelDisc(;potential::MestelPotential=MestelPotential(),q::IntorFloat=11.44,G::Float64=1.0)
    return MestelDisc(potential,q,G)
end

"""
    MestelDistribution(EL::Tuple{Float64,Float64},df::MestelDisc)
Mestel distribution function.
"""
function DistributionFunction(EL::Tuple{Float64,Float64},df::MestelDisc)::Float64

    return MestelDistribution(EL,df)
end

"""
    dFdE(EL::Tuple{Float64,Float64},df::MestelDisc)
Mestel DF derivative w.r.t. E.
"""
function DFDE(EL::Tuple{Float64,Float64},df::MestelDisc)::Float64
    return MesteldFdE(EL,df)
end

"""
    dFdL(E, L[, C, q, sigma])
Mestel DF derivative w.r.t. E.
"""
function DFDL(EL::Tuple{Float64,Float64},df::MestelDisc)::Float64
    return MesteldFdL(EL,df)
end

