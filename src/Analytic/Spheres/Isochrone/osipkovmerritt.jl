


"""
OsipkovMerrittIsochrone([potential])

Osipkov-Merritt anisotropy radius Isochrone distribution function. Uses OrbitalElements.NumericalIsochrone by default.
"""
function OsipkovMerrittIsochroneEL(ra::Float64; potential::IsochronePotential=NumericalIsochrone())
    return OsipkovMerrittIsochroneEL(ra,potential)
end
function OsipkovMerrittIsochroneJL(ra::Float64; potential::IsochronePotential=NumericalIsochrone())
    return OsipkovMerrittIsochroneJL(ra,potential)
end

"""
    osipkovmerritt_Q(EL::Tuple{Float64,Float64}, df::OsipkovMerrittIsochrone)

Q(E,L), the variable for anisotropy distribution functions.

when ra->infty, this is isotropic.
"""
function osipkovmerritt_Q(EL::Tuple{Float64,Float64}, df::OsipkovMerrittIsochrone)
    E,L = EL

    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc

    Q = (E + (L^(2)*(df.potential.bc^2))/(2.0*df.ra^(2)))/scaleEnergy
    return Q # Output

end


"""
jacobian for converting dF/dQ dQ/dE -> dF/dE
"""
function osipkovmerritt_dQdE(EL::Tuple{Float64,Float64}, df::OsipkovMerrittIsochrone)
    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    return 1.0/scaleEnergy # Output
end

"""
jacobian for converting dF/dQ dQ/dL -> dF/dL
"""
function osipkovmerritt_dQdL(EL::Tuple{Float64,Float64}, df::OsipkovMerrittIsochrone)
    E,L = EL
    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    return (L*(df.potential.bc^2))/(scaleEnergy*df.ra^(2)) # Output
end


"""
Saha distribution function
ra is the anisotropy radius
"""
function DistributionFunction(EL::Tuple{Float64,Float64}, df::OsipkovMerrittIsochrone)::Float64

    Q       = osipkovmerritt_Q(EL, df)
    scaleDF = dfscale(df)

    # dimensionless anisotropy radius
    gamma = (df.potential.bc/df.ra)^2

    prefactor = (df.potential.M/scaleDF) * (1/(128*sqrt(2)*(pi^3))) * (1/((1-Q)^4))

    term1 = sqrt(Q)*(27 + 77gamma - (66+286gamma)*Q + (320+136gamma)*Q*Q -(240+32gamma)*Q*Q*Q + 64*Q*Q*Q*Q)
    term2 = ((3*asin(sqrt(Q)))/sqrt((1-Q))) * ((-9 + 17gamma) + (28-44gamma)*Q + (16-8gamma)*Q*Q)

    return prefactor * (term1+term2)
end


"""
Saha distribution function derivative w.r.t. Q
ra is the anisotropy radius

"""
function osipkovmerritt_dFdQ(EL::Tuple{Float64,Float64}, df::OsipkovMerrittIsochrone)
    
    Q = osipkovmerritt_Q(EL,df)
    scaleDF = dfscale(df)

    # dimensionless anisotropy radius
    gamma = (df.potential.bc/df.ra)^2

    prefactor = df.potential.M/scaleDF

    prefactor2 = 1/(128*sqrt(2)*(pi^3))
    prefactor3 = - (1/(2Q*(1-Q)^(11/2)))

    term1 = (15Q*(5-13gamma - 52Q+ 68gamma*Q + 8*(-2+gamma)*(Q^2)))*asin(sqrt(Q))
    term2 = sqrt(Q) * (sqrt(1-Q)*(-27 + gamma * (-77 + Q*(319+2Q*(375+4*Q*(-23 + 4*Q)))) + Q*(9+2*Q*(-635+8*Q*(45+Q*(-21+4*Q)))) ) - 3 * (-1+Q)*(9-4*Q*(7+4*Q)+gamma*(-17+44*Q+8*(Q^2))*(1/(1-Q))))

    return prefactor*prefactor2*prefactor3*(term1+term2)

end

function DFDE(EL::Tuple{Float64,Float64}, df::OsipkovMerrittIsochrone)::Float64

    Q = osipkovmerritt_Q(EL,df)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTEN TION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    # Value of dF/dQ
    dFdQ = osipkovmerritt_dFdQ(EL,df)

    # Value of dQ/dE
    dQdE = osipkovmerritt_dQdE(EL,df)

    return dFdQ * dQdE
end

function DFDL(EL::Tuple{Float64,Float64}, df::OsipkovMerrittIsochrone)::Float64

    Q = osipkovmerritt_Q(EL,df)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTEN TION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    # Value of dF/dQ
    dFdQ = osipkovmerritt_dFdQ(EL,df)

    # Value dQ/dL
    dQdL = osipkovmerritt_dQdL(EL,df)

    return dFdQ * dQdL
end
