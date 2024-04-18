
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
function F(EL::Tuple{Float64,Float64}, df::IsotropicPlummer)
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
function dFdE(E::Float64,df::IsotropicPlummer)

    scaleEnergy = - df.potential.G * df.potential.M / df.potential.bc
    scaleDF     = dfscale(df)

    prefactor = df.potential.M * scaleDF * (12*sqrt(2)/(pi^3))

    return - prefactor * (E/scaleEnergy)^(5/2)

end


"""
calculate n.dF/dE, for the isotropic plummer case.
"""
function ndFdJ(EL::Tuple{Float64,Float64},ΩΩ::Tuple{Float64,Float64},resonance::Resonance, df::IsotropicPlummer)
    E,L = EL
    Ω1,Ω2 = ΩΩ
    ndotΩ = resonance.number[1]*Ω1 + resonance.number[2]*Ω2

    dFdE = dFdE(E,df)

    # Value of n.dF/dJ
    result = ndotOmega*dFdE

    return result
end