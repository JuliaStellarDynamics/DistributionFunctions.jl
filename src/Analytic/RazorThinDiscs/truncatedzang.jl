#####
# Truncated full tapered DF and derivatives
#
# Adding the truncation: no particles beyond Rmax
#####

"""
TruncatedZangDisc([potential])

Zang disc distribution function.
"""
function TruncatedZangDisc(;potential::MestelPotential=MestelPotential(),q::IntorFloat=11.44,ν::Int64=4,Rin::Float64=1.0,μ::Int64=5,Rout::Float64=11.5,Rmax::Float64=20.,G::Float64=1.0)
    return TruncatedZangDisc(potential,q,ν,Rin,μ,Rout,Rmax,G)
end


function truncationcriteria(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)::Float64

    E,L = EL

    if (L <= 0.) || (L > df.Rmax*df.potential.V0) || (E < (df.potential.V0^2)/2 + ψ(L/df.potential.V0,df.potential)) || (E > ψ(df.Rmax,df.potential) + L^2/(2*df.Rmax^2))
        return 0.
    end

    # if not returning zero, return one
    return 1.0
end

"""
    TruncatedZangDistribution(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)

Zang star distribution function enforcing ra <= Rmax.
"""
function Distribution(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)::Float64

    # will return zero if outside truncated region
    result = truncationcriteria(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)

    return result * MestelDistribution(EL,df) * ZangOuterTaper(EL,df) * ZangInnerTaper(EL,df)
end


"""
    dFdE(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)
Zang star DF derivative w.r.t. E.
"""
function DFDE(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)::Float64
    # will return zero if outside truncated region
    result = truncationcriteria(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)

    return result * MesteldFdE(EL,df) * ZangOuterTaper(EL,df) * ZangInnerTaper(EL,df)
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)
Zang star DF derivative w.r.t. L.
"""
function DFDL(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)::Float64

    mesDF = MestelDistribution(EL,df)
    intap = ZangInnerTaper(EL,df)
    outap = ZangOuterTaper(EL,df)

    # will return zero if outside truncated region
    result = truncationcriteria(EL::Tuple{Float64,Float64},df::TruncatedZangDisc)

    return result * ( MesteldFdL(EL,df) * intap * outap + 
            mesDF * ZangInnerTaperdL(EL,df) * outap + 
            mesDF * intap * ZangOuterTaperdL(EL,df) )
end
