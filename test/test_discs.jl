

# define our resonance of interest for mapping
n1,n2 = 1,2
a,e = 0.1,0.0

@testset "discdistributionfunctions" begin
    @testset "mestel" begin
        DDF = MestelDisc()
        EL = EL_from_ae(a,e,DDF.potential)
        ΩΩ = frequencies_from_ae(a,e,DDF.potential)
        resonance = Resonance(n1,n2,DDF.potential)
        @test DistributionFunction(EL, DDF) ≈ 4.397057 atol=1e-6
        @test gradient(EL, DDF)[1] ≈ -54.699389 atol=1e-6
        @test gradient(EL, DDF)[2] ≈ 503.0233217077802 atol=1e-6
    end
    @testset "zang" begin
        DDF = ZangDisc()
        EL = EL_from_ae(a,e,DDF.potential)
        ΩΩ = frequencies_from_ae(a,e,DDF.potential)
        resonance = Resonance(n1,n2,DDF.potential)
        @test DistributionFunction(EL, DDF) ≈ 0.0004396617345 atol=1e-6
        @test gradient(EL, DDF)[1] ≈ -5.468845093620248e-7 atol=1e-6
        @test gradient(EL, DDF)[2] ≈ 0.067882013 atol=1e-6
    end
    @testset "truncatedzang" begin
        # test inside the valid region
        DDF = TruncatedZangDisc()
        DFcomp = ZangDisc()
        EL = EL_from_ae(a,e,DDF.potential)
        ΩΩ = frequencies_from_ae(a,e,DDF.potential)
        resonance = Resonance(n1,n2,DDF.potential)
        @test DistributionFunction(EL, DDF) == DistributionFunction(EL, DFcomp)
        @test gradient(EL, DDF)[1] == gradient(EL, DFcomp)[1]
        @test gradient(EL, DDF)[2] == gradient(EL, DFcomp)[2]
        # test outside the valid region?
    end
end