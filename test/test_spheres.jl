# Does not test :
#   - 
# @IMPROVE add undefined scale functions tests once defined

using OrbitalElements

# define our resonance of interest for mapping
n1,n2 = -1,2

@testset "distributionfunctions" begin
    @testset "isochrone" begin
        @testset "isotropic" begin
            IsoDF = IsotropicIsochrone()
            EL = EL_from_ae(1.0,0.5,IsoDF.potential)
            resonance = Resonance(n1,n2,IsoDF.potential)
            @test F(EL, IsoDF) ≈ 0.023905 atol=1e-6
        end
        @testset "osipkovmerritt" begin
            
        end
    end
    @testset "plummer" begin
    end
end