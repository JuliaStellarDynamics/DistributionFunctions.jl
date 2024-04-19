# Does not test :
#   - 
# @IMPROVE add undefined scale functions tests once defined

using OrbitalElements

# define our resonance of interest for mapping
n1,n2 = -1,2

@testset "distributionfunctions" begin
    @testset "isochrone" begin
        IsoDF = IsotropicIsochrone()
        EL = EL_from_ae(1.0,0.5,IsoDF.potential)
        resonance = Resonance(n1,n2,IsoDF.potential)
    end
end