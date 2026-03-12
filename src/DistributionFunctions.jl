module DistributionFunctions

#####################################
# Dependencies
#####################################
using OrbitalElements             # potentials


#####################################
# Exports
#####################################
export DistributionFunction

# functions common to all DistributionFunction
export gradient

# types for multiple dispatch
export ErgodicDF,EnergyAngularMomentumDF,ActionDF

# spheres
export PlummerDF,IsotropicPlummer,OsipkovMerrittPlummer,OsipkovMerrittPlummerEL,OsipkovMerrittPlummerJL
export IsochroneDF,IsotropicIsochrone,OsipkovMerrittIsochroneEL,OsipkovMerrittIsochroneJL

# discs
export MestelDisc,ZangDisc,TruncatedZangDisc
export ToomreDisc, ToomreDiscOdd


#####################################
# Includes
#####################################
include("Generic/DFs.jl")

end # module