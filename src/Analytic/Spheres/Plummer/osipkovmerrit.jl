

"""
OsipkovMerrittPlummer([potential])

Osipkov-Merritt anisotropy radius Plummer distribution function. Uses OrbitalElements.NumericalPlummer by default.
"""
function OsipkovMerrittPlummerEL(ra::Float64; potential::PlummerPotential=NumericalPlummer())
    return OsipkovMerrittPlummerEL(ra,potential)
end
function OsipkovMerrittPlummerJL(ra::Float64; potential::PlummerPotential=NumericalPlummer())
    # this doesn't exist. use (E,L) instead
    return OsipkovMerrittPlummerJL(ra,potential)
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
function DistributionFunction(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)

    Q       = osipkovmerritt_Q(EL, df)
    scaleDF = dfscale(df)

    # dimensionless anisotropy radius
    gamma = (df.potential.bc/df.ra)^2

    prefactor = (df.potential.M/scaleDF) * (24*sqrt(2)/(7*(pi^3)))

    return prefactor * (Q^(7/2)) * (1-gamma+7gamma/(16*(Q^2)))

end

function DFDE(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)::Float64

    Q = osipkovmerritt_Q(EL, df)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    # Value of dF/dQ
    dFdQ = osipkovmerritt_dFdQ(EL,df)

    # Value of dQ/dE, dQ/dL
    dQdE = osipkovmerritt_dQdE(EL,df)

    return dFdQ*dQdE
end

function DFDL(EL::Tuple{Float64,Float64}, df::OsipkovMerrittPlummer)::Float64

    Q = osipkovmerritt_Q(EL, df)

    # If Q is outside of the [0,1]--range, we set the function to 0.0
    # ATTENTION, this is a lazy implementation -- it would have been much better to restrict the integration domain
    if (!(0.0 <= Q <= 1.0)) # If Q is outside of the [0,1]-range, we set the function to 0
        return 0.0 # Outside of the physically allowed orbital domain
    end

    # Value of dF/dQ
    dFdQ = osipkovmerritt_dFdQ(EL,df)

    # Value of  dQ/dL
    dQdL =  osipkovmerritt_dQdL(EL,df)

    return dFdQ*dQdL
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

