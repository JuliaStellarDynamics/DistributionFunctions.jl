

"""
MestelDisc([potential])

Mestel disc distribution function.
"""
function MestelDisc(;potential::MestelPotential=MestelPotential(),q::IntorFloat=11.44,ξDF::Float64=1.0,G::Float64=1.0)
    return MestelDisc(potential,q,ξDF,G)
end

"""
    MestelDF(EL::Tuple{Float64,Float64},df::MestelDisc)
Mestel distribution function.
"""
function F(EL::Tuple{Float64,Float64},df::MestelDisc)::Float64

    return MestelDF(EL,df)
end

"""
    dFdE(EL::Tuple{Float64,Float64},df::MestelDisc)
Mestel DF derivative w.r.t. E.
"""
function dFdE(EL::Tuple{Float64,Float64},df::MestelDisc)::Float64
    return MesteldFdE(EL,df)
end

"""
    dFdL(E, L[, C, q, sigma])
Mestel DF derivative w.r.t. E.
"""
function dFdL(EL::Tuple{Float64,Float64},df::MestelDisc)::Float64
    return MesteldFdL(EL,df)
end


"""
    ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::MestelDisc)
Zang star DF derivative w.r.t. the actions J.
"""
function ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::MestelDisc)::Float64

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