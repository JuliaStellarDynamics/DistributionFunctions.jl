
"""
IsotropicPlummer([potential])

Isotropic Plummer distribution function. Uses OrbitalElements.NumericalPlummer by default.
"""
function IsotropicPlummer(;potential::PlummerPotential=NumericalPlummer())
    return IsotropicPlummer(potential)
end



"""
The distribution function for the isotropic Plummer model
"""
function Distribution(EL::Tuple{Float64,Float64}, df::IsotropicPlummer)
    E,L = EL

    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    scaleDF     = dfscale(df)

    prefactor = (df.potential.M/scaleDF) * (24*sqrt(2)/(7*(pi^3)))

    return prefactor * (E/scaleEnergy)^(7/2)

end

"""
the derivative of the isotropic distribution function

F[EN_] = (24 Sqrt[2])/(7 Pi^3) * EN^(7/2)
D[F[EN], {EN}]
"""
function DFDE(EL::Tuple{Float64,Float64},df::IsotropicPlummer)

    E,L = EL

    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    scaleDF     = dfscale(df)

    prefactor = df.potential.M * scaleDF * (12*sqrt(2)/(pi^3))

    return - prefactor * (E/scaleEnergy)^(5/2)

end

function DFDL(EL::Tuple{Float64,Float64},df::IsotropicPlummer)::Float64

    return 0.0

end

