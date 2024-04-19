
"""
IsotropicIsochrone([potential])

Isotropic Isochrone distribution function. Uses OrbitalElements.NumericalIsochrone by default.
"""
function IsotropicIsochrone(;potential::IsochronePotential=NumericalIsochrone())
    return IsotropicIsochrone(potential)
end



"""isochroneisotropicDF(E[,bc,M,G])
the isotropic DF
"""
function F(EL::Tuple{Float64,Float64}, df::IsotropicIsochrone)
    E,L = EL

    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    scaleDF     = dfscale(df)

    # rescale dimensionless energy to goes between 0 and 1/2
    mE = E/scaleEnergy

    return (sqrt(mE)*df.potential.M*(27.0+2.0*mE*(-1.0+4.0*mE)*(33.0+4.0*mE*(-7.0+2.0*mE))+
    (3.0*(-9.0+4.0*mE*(7.0+4.0*mE))*asin(sqrt(mE)))/sqrt(-((-1.0+mE)*mE))))/
    (128.0*sqrt(2.0)*(-1.0+mE)^(4)*scaleDF*(pi)^(3))

end


"""isochroneisotropicdDFdE(E[,bc,M,G])
the isotropic DF derivative w.r.t. E
"""
function dFdE(EL::Tuple{Float64,Float64}, df::IsotropicIsochrone)

    E,L = EL

    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    scaleDF     = dfscale(df)

    # rescale dimensionless energy to goes between 0 and 1/2
    mE = E/scaleEnergy

    # return the isotropic derivative
    return (df.potential.M*((-1.0+mE)*sqrt(mE)*(-75.0+2.0*mE*(-659.0+8.0*mE*(45.0+mE*(-21.0+4.0*mE))))+
    15.0*sqrt(1.0-mE)*(-5.0+4.0*mE*(13.0+4.0*mE))*asin(sqrt(mE))))/
    (256.0*sqrt(2.0)*scaleEnergy*(-1.0+mE)^(6)*scaleDF*(pi)^(3))

end


function ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::IsotropicIsochrone)
    Ω1,Ω2 = ΩΩ
    ndotΩ = resonance.number[1]*Ω1 + resonance.number[2]*Ω2

    dFdEval = dFdE(EL,df)

    # Value of n.dF/dJ
    result = ndotΩ*dFdEval

    return result
end
