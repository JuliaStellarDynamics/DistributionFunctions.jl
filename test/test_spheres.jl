# Does not test :
#   - 
# @IMPROVE add undefined scale functions tests once defined

# define our resonance of interest for mapping
n1,n2 = -1,2
a,e = 1.0,0.5

@testset "sphericaldistributionfunctions" begin
    @testset "isochrone" begin
        @testset "isotropic" begin
            IsoDF = IsotropicIsochrone()
            EL = EL_from_ae(a,e,IsoDF.potential)
            ΩΩ = frequencies_from_ae(a,e,IsoDF.potential)
            resonance = Resonance(n1,n2,IsoDF.potential)
            @test DistributionFunction(EL, IsoDF) ≈ 0.023905 atol=1e-6
            @test gradient(EL, IsoDF)[1] ≈ -0.315496 atol=1e-6
            @test gradient(EL, IsoDF)[2] == 0.0
        end
        @testset "osipkovmerrittEL" begin
            ra = 1.0
            OMDF = OsipkovMerrittIsochroneEL(ra)
            EL = EL_from_ae(1.0,0.5,OMDF.potential)
            ΩΩ = frequencies_from_ae(a,e,OMDF.potential)
            resonance = Resonance(n1,n2,OMDF.potential)
            @test DistributionFunction(EL, OMDF) ≈ 0.020028 atol=1e-6
            @test gradient(EL, OMDF)[1] ≈ -0.07563618 atol=1e-6
            @test gradient(EL, OMDF)[2] ≈ -0.01926609 atol=1e-6        
        end
    end
    @testset "plummer" begin
        @testset "isotropic" begin
            IsoDF = IsotropicPlummer()
            EL = EL_from_ae(1.0,0.5,IsoDF.potential)
            ΩΩ = frequencies_from_ae(a,e,IsoDF.potential)
            resonance = Resonance(n1,n2,IsoDF.potential)
            @test DistributionFunction(EL, IsoDF) ≈ 0.015042 atol=1e-6
            @test gradient(EL, IsoDF)[1] ≈ -0.1027822 atol=1e-6
            @test gradient(EL, IsoDF)[2] == 0.0
        end
        @testset "osipkovmerritt" begin
            ra = 0.75
            OMDF = OsipkovMerrittPlummerEL(ra)
            EL = EL_from_ae(1.0,0.5,OMDF.potential)
            ΩΩ = frequencies_from_ae(a,e,OMDF.potential)
            resonance = Resonance(n1,n2,OMDF.potential)
            @test DistributionFunction(EL, OMDF) ≈ 0.021509 atol=1e-6  
            @test gradient(EL, OMDF)[1] ≈ -0.07755405 atol=1e-6
            @test gradient(EL, OMDF)[2] ≈ -0.060270980 atol=1e-6
        end
    end
end