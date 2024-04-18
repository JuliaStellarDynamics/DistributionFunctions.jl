module DistributionFunctions

#####################################
# Dependencies
#####################################
using OrbitalElements             # potentials, resonances


#####################################
# Exports
#####################################
export DistributionFunction
export PlummerDistributionFunction,IsotropicPlummer,OsipkovMerrittPlummer

export F,ndFdJ
export dfscale

#####################################
# Includes
#####################################
include("Generic/DFs.jl")

end # module