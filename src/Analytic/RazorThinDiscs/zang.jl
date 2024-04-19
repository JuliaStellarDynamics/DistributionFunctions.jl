
"""
ZangDisc([potential])

Zang disc distribution function.
"""
function ZangDisc(;potential::MestelPotential=MestelPotential(),q::IntorFloat=11.44,ν::Int64=4,Rin::Float64=1.0,μ::Int64=5,Rout::Float64=11.5,ξDF::Float64=1.0,G::Float64=1.0)
    return ZangDisc(potential,q,ν,Rin,μ,Rout,ξDF,G)
end


#####
# Zang tapering and derivatives
#####

"""
    ZangInnerTaper(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang inner tapering.
"""
function ZangInnerTaper(EL::Tuple{Float64,Float64},df::ZangDistributionFunction)::Float64

    E,L = EL
    return (L^df.ν) / ( (df.Rin*df.potential.V0)^(df.ν) + (L^df.ν) )
end

"""
    ZangInnerTaperdL(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang inner tapering derivative.
"""
function ZangInnerTaperdL(EL::Tuple{Float64,Float64},df::ZangDistributionFunction)::Float64

    E,L = EL
    return (df.Rin*df.potential.V0)^(df.ν) * df.ν * (L)^(df.ν-1) / ( (df.Rin*df.potential.V0)^(df.ν) + (L)^(df.ν) )^(2)
end

"""
    Zang_outer_tapering(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang outer tapering.
"""
function ZangOuterTaper(EL::Tuple{Float64,Float64},df::ZangDistributionFunction)::Float64

    E,L = EL
    return (df.Rout*df.potential.V0)^(df.μ) / ( (df.Rout*df.potential.V0)^(df.μ) + (L^df.μ) )
end
"""
    Zang_outer_tapering_dL(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang outer tapering derivative.
"""
function ZangOuterTaperdL(EL::Tuple{Float64,Float64},df::ZangDistributionFunction)::Float64

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
function F(EL::Tuple{Float64,Float64},df::ZangDisc)::Float64
    return MestelDF(EL,df) * ZangOuterTaper(EL,df) * ZangInnerTaper(EL,df)
end

"""
    dFdE(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang star DF derivative w.r.t. E.
"""
function dFdE(EL::Tuple{Float64,Float64},df::ZangDistributionFunction)::Float64
    return MesteldFdE(EL,df) * ZangOuterTaper(EL,df) * ZangInnerTaper(EL,df)
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::ZangDisc)
Zang star DF derivative w.r.t. L.
"""
function dFdL(EL::Tuple{Float64,Float64},df::ZangDistributionFunction)::Float64

    mesDF = MestelDF(EL,df)
    intap = ZangInnerTaper(EL,df)
    outap = ZangOuterTaper(EL,df)
    return ( MesteldFdL(EL,df) * intap * outap + 
            mesDF * ZangInnerTaperdL(EL,df) * outap + 
            mesDF * intap * ZangOuterTaperdL(EL,df) )
end

"""
    ZangndDFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::ZangDisc)
Zang star DF derivative w.r.t. the actions J.
"""
function ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::ZangDisc)::Float64

    E,L = EL
    Ω1,Ω2 = ΩΩ
    n1,n2 = resonance.number
    ndotΩ = n1*Ω1 + n2*Ω2

    if L <= 0.
        println("WARNING: L <= 0.")
        return 0.
    end

    dDFdE = dFdE(EL,df)
    dDFdL = dFdL(EL,df)
    
    return df.ξDF * (ndotΩ*dDFdE + n2*dDFdL)
end

