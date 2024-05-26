module DistributionFunctions

#####################################
# Dependencies
#####################################
using OrbitalElements             # potentials, resonances


#####################################
# Exports
#####################################
export DistributionFunction
export Distribution,DFDE,DFDL,gradient

# types for multiple dispatch
export EnergyOnlyDistributionFunction,EnergyAngularMomentumDistributionFunction,ActionDistributionFunction

# spheres
export PlummerDistributionFunction,IsotropicPlummer,OsipkovMerrittPlummerEL,OsipkovMerrittPlummerJL
export IsochroneDistributionFunction,IsotropicIsochrone,OsipkovMerrittIsochroneEL,OsipkovMerrittIsochroneJL

# discs
export MestelDisc,ZangDisc,TruncatedZangDisc


#####################################
# Includes
#####################################
include("Generic/DFs.jl")

end # module