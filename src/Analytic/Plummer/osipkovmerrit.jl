

"""
OsipkovMerrittPlummer([potential])

Osipkov-Merritt anisotropy radius Plummer distribution function. Uses OrbitalElements.NumericalPlummer by default.
"""
function OsipkovMerrittPlummer(ra::Float64; potential::PlummerPotential=NumericalPlummer())
    return OsipkovMerrittPlummer(ra,potential)
end


"""
    osipkovmerritt_Q(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)

Q(E,L), the variable for anisotropy distribution functions.

when ra->infty, this is isotropic.
"""
function osipkovmerritt_Q(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)
    E,L = EL

    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc

    Q = (E + (L^(2)*(df.potential.bc^2))/(2.0*df.ra^(2)))/scaleEnergy
    return Q # Output

end

"""
jacobian for converting dF/dQ dQ/dE -> dF/dE
"""
function osipkovmerritt_dQdE(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)
    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    return 1.0/scaleEnergy # Output
end

"""
jacobian for converting dF/dQ dQ/dL -> dF/dL
"""
function osipkovmerritt_dQdL(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)
    E,L = EL
    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    return (L*(df.potential.bc^2))/(scaleEnergy*df.ra^(2)) # Output
end


"""
the anisotropic distribution function from PlummerPlus
"""
function F(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)

    Q       = osipkovmerritt_Q(EL, df)
    scaleDF = dfscale(df)

    # dimensionless anisotropy radius
    gamma = (df.potential.bc/df.ra)^2

    prefactor = (df.potential.M/scaleDF) * (24*sqrt(2)/(7*(pi^3)))

    return prefactor * (Q^(7/2)) * (1-gamma+7gamma/(16*(Q^2)))

end


"""
the derivative of the anisotropic distribution function
"""
function osipkovmerritt_dFdQ(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)

    Q       = osipkovmerritt_Q(EL, df)
    scaleDF = dfscale(df)

    # dimensionless anisotropy radius: ALREADY SQUARED
    gamma = (df.potential.bc/df.ra)^2

    prefactor = (df.potential.M/scaleDF) * ((3*sqrt(2))/(4*(pi^3)))
    return prefactor * sqrt(Q) * ( (gamma)*(3 - 16(Q^2)) + 16(Q^2)) 

end


"""
calculate n.dF/dJ, for the osipkov-merritt plummer case.

@IMPROVE, make a distribution function fix. perhaps only integrate in regions where the DF is positive?
"""
function ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::OsipkovMerrittPlummer)

    Ω1,Ω2 = ΩΩ
    n1,n2 = resonance.number[1],resonance.number[2]
    ndotΩ = n1*Ω1 + n2*Ω2

    Q = osipkovmerritt_Q(EL, df)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    # Value of dF/dQ
    dFdQ = osipkovmerritt_dFdQ(EL,df)

    # Values of dQ/dE, dQ/dL
    dQdE, dQdL = osipkovmerritt_dQdE(EL,df), osipkovmerritt_dQdL(EL,df)

    # Value of n.dF/dJ
    result = dFdQ*(dQdE*ndotΩ + n2*dQdL)

    return result

end
