
"""
ZangDisc([potential])

Zang disc distribution function.
"""
function ZangDisc(;potential::MestelPotential=MestelPotential(),q::IntorFloat=11.44,ν::Int64=4,Rin::Float64=1.0,μ::Int64=5,Rout::Float64=11.5,G::Float64=1.0)
    return ZangDisc(potential,q,ν,Rin,μ,Rout,G)
end


#####
# Zang tapering and derivatives
#####

"""
    ZangInnerTaper(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang inner tapering.
"""
function ZangInnerTaper(EL::Tuple{Float64,Float64},df::ZangDF)::Float64

    E,L = EL
    return (L^df.ν) / ( (df.Rin*df.potential.V0)^(df.ν) + (L^df.ν) )
end

"""
    ZangInnerTaperdL(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang inner tapering derivative.
"""
function ZangInnerTaperdL(EL::Tuple{Float64,Float64},df::ZangDF)::Float64

    E,L = EL
    return (df.Rin*df.potential.V0)^(df.ν) * df.ν * (L)^(df.ν-1) / ( (df.Rin*df.potential.V0)^(df.ν) + (L)^(df.ν) )^(2)
end

"""
    Zang_outer_tapering(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang outer tapering.
"""
function ZangOuterTaper(EL::Tuple{Float64,Float64},df::ZangDF)::Float64

    E,L = EL
    return (df.Rout*df.potential.V0)^(df.μ) / ( (df.Rout*df.potential.V0)^(df.μ) + (L^df.μ) )
end
"""
    Zang_outer_tapering_dL(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang outer tapering derivative.
"""
function ZangOuterTaperdL(EL::Tuple{Float64,Float64},df::ZangDF)::Float64

    E,L = EL
    return - (df.Rout*df.potential.V0)^(df.μ) * df.μ * (L)^(df.μ-1) / ( (df.Rout*df.potential.V0)^(df.μ) + (L)^(df.μ) )^(2)
end

#####
# Full tapered DF and derivatives
#####

"""
    F(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang star distribution function.
"""
function DistributionFunction(EL::Tuple{Float64,Float64},df::ZangDisc)::Float64
    return MestelDistribution(EL,df) * ZangOuterTaper(EL,df) * ZangInnerTaper(EL,df)
end

"""
    dFdE(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang star DF derivative w.r.t. E.
"""
function DFDE(EL::Tuple{Float64,Float64},df::ZangDF)::Float64
    return MesteldFdE(EL,df) * ZangOuterTaper(EL,df) * ZangInnerTaper(EL,df)
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang star DF derivative w.r.t. L.
"""
function DFDL(EL::Tuple{Float64,Float64},df::ZangDF)::Float64

    mesDF = MestelDistribution(EL,df)
    intap = ZangInnerTaper(EL,df)
    outap = ZangOuterTaper(EL,df)
    return ( MesteldFdL(EL,df) * intap * outap + 
            mesDF * ZangInnerTaperdL(EL,df) * outap + 
            mesDF * intap * ZangOuterTaperdL(EL,df) )
end
