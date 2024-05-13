

# define our resonance of interest for mapping
n1,n2 = 1,2
a,e = 0.1,0.0

@testset "discdistributionfunctions" begin
    @testset "mestel" begin
        DDF = MestelDisc()
        EL = EL_from_ae(a,e,DDF.potential)
        ΩΩ = frequencies_from_ae(a,e,DDF.potential)
        resonance = Resonance(n1,n2,DDF.potential)
        @test Distribution(EL, DDF) ≈ 4.397057 atol=1e-6
        @test DFDE(EL, DDF) ≈ -54.699389 atol=1e-6
        @test DFDL(EL, DDF) ≈ 503.0233217077802 atol=1e-6
        @test ndFdJ(EL,ΩΩ,resonance, DDF) ≈ -861.50732 atol=1e-6
    end
    @testset "zang" begin
        DDF = ZangDisc()
        EL = EL_from_ae(a,e,DDF.potential)
        ΩΩ = frequencies_from_ae(a,e,DDF.potential)
        resonance = Resonance(n1,n2,DDF.potential)
        @test Distribution(EL, DDF) ≈ 0.0004396617345 atol=1e-6
        @test DFDE(EL, DDF) ≈ -5.468845093620248e-7 atol=1e-6
        @test DFDL(EL, DDF) ≈ 0.067882013 atol=1e-6
        @test ndFdJ(EL,ΩΩ,resonance, DDF) ≈ 0.135745354 atol=1e-6
    end
    @testset "truncatedzang" begin
        # test inside the valid region
        DDF = TruncatedZangDisc()
        DFcomp = ZangDisc()
        EL = EL_from_ae(a,e,DDF.potential)
        ΩΩ = frequencies_from_ae(a,e,DDF.potential)
        resonance = Resonance(n1,n2,DDF.potential)
        @test Distribution(EL, DDF) == Distribution(EL, DFcomp)
        @test DFDE(EL, DDF) == DFDE(EL, DFcomp)
        @test DFDL(EL, DDF) == DFDL(EL, DFcomp)
        @test ndFdJ(EL,ΩΩ,resonance, DDF) == ndFdJ(EL,ΩΩ,resonance, DFcomp)
        # test outside the valid region?
    end
end