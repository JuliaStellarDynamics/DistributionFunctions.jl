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
export distribution,gradient

# types for multiple dispatch
export ErgodicDF,EnergyAngularMomentumDF,ActionDF

# spheres
# change to DF
export PlummerDF,IsotropicPlummer,OsipkovMerrittPlummerEL,OsipkovMerrittPlummerJL
export IsochroneDF,IsotropicIsochrone,OsipkovMerrittIsochroneEL,OsipkovMerrittIsochroneJL

# discs
# add DF here
export MestelDisc,ZangDisc,TruncatedZangDisc


#####################################
# Includes
#####################################
include("Generic/DFs.jl")

end # module