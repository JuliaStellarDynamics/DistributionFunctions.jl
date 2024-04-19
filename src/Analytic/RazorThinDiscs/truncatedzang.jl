#####
# Truncated full tapered DF and derivatives
#
# Adding the truncation: no particles beyond Rmax
#####

"""
TruncatedZangDisc([potential])

Zang disc distribution function.
"""
function TruncatedZangDisc(;potential::MestelPotential=MestelPotential(),q::IntorFloat=11.44,ν::Int64=4,Rin::Float64=1.0,μ::Int64=5,Rout::Float64=11.5,Rmax::Float64=20.,ξDF::Float64=1.0,G::Float64=1.0)
    return TruncatedZangDisc(potential,q,ν,Rin,μ,Rout,Rmax,ξDF,G)
end



"""
    TruncatedZangDF(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)

Zang star distribution function enforcing ra <= Rmax.
"""
function F(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)::Float64

    E,L = EL

    if (L <= 0.) || (L > df.Rmax*df.potential.V0) || (E < (df.potential.V0^2)/2 + ψ(L/df.potential.V0,df.potential)) || (E > ψ(df.Rmax,df.potential) + L^2/(2*df.Rmax^2))
        return 0.
    end
    return MestelDF(EL,df) * ZangOuterTaper(EL,df) * ZangInnerTaper(EL,df)
end



"""
    TruncatedZangndDFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::TruncatedZangDisc)

Zang star DF derivative w.r.t. the actions J enforcing ra <= Rmax.
"""
function ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::TruncatedZangDisc)::Float64

    E,L = EL
    Ω1,Ω2 = ΩΩ
    n1,n2 = resonance.number
    ndotΩ = n1*Ω1 + n2*Ω2

    if (L <= 0.) || (L > df.Rmax*df.potential.V0) || (E < (df.potential.V0^2)/2 + ψ(L/df.potential.V0,df.potential)) || (E > ψ(df.Rmax,df.potential) + L^2/(2*df.Rmax^2))
        return 0.
    end

    dDFdE = dFdE(EL,df)
    dDFdL = dFdL(EL,df)
    
    return df.ξDF * (ndotΩ*dDFdE + n2*dDFdL)
end