module DistributionFunctions

#####################################
# Dependencies
#####################################
using OrbitalElements             # potentials, resonances


#####################################
# Exports
#####################################
export DistributionFunction
export F,ndFdJ

# spheres
export PlummerDistributionFunction,IsotropicPlummer,OsipkovMerrittPlummer
export IsochroneDistributionFunction,IsotropicIsochrone,OsipkovMerrittIsochrone
export dfscale

# discs
export MestelDisc,ZangDisc,TruncatedZangDisc
export σMestelDF,NormConstMestelDF,MestelDF,MesteldFdE,MesteldFdL


#####################################
# Includes
#####################################
include("Generic/DFs.jl")

end # module