using OrbitalElements

abstract type DistributionFunction end

abstract type PlummerDistributionFunction <: DistributionFunction end

struct IsotropicPlummer <: PlummerDistributionFunction
    potential::PlummerPotential # Potential model
end

function IsotropicPlummer(;potential::PlummerPotential=NumericalPlummer())
    return IsotropicPlummer(potential)
end

"""
the Plummer distribution function scale
"""
function dfscale(df::PlummerDistributionFunction)
    return (df.potential.G*df.potential.M*df.potential.bc)^(-3/2)
end

IPlum = IsotropicPlummer()
IPlum = IsotropicPlummer(SemiAnalyticPlummer())

EL = EL_from_ae(1.0,0.5,IPlum.potential)

F(EL, IPlum)

APlum = OsipkovMerrittPlummer(1.0)


n1,n2 = -1,2
resonance = Resonance(n1,n2,IPlum.potential)