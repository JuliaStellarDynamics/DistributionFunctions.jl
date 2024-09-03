
# a, c, G, M in potential

struct KuzminSpheroid{modelT<:KuzminPotentials} <: KuzminPotentialDF
    potential::modelT # potential model
end




"""
   KuzminDistribution(ELz::Tuple{Float64,Float64},df::MestelDisc)
Kuzmin-Kutuzov ddistribution function.
"""
function KuzminDistribution(ELz::Tuple{Float64,Float64}, df::KuzminSpheroid)::Float64

    E,Lz = ELz
    a = df.potential.a
    c = df.potential.c
    G = df.potential.G
    M = df.potential.M
    
    # todo
    
end

function DFDE(ELz::Tuple{Float64,Float64}, df::KuzminSpheroid)::Float64

    E,Lz = ELz
    a = df.potential.a
    c = df.potential.c
    G = df.potential.G
    M = df.potential.M

    # todo
    
end

function DFDLz(ELz::Tuple{Float64,Float64}, df::KuzminSpheroid)::Float64

    E,Lz = ELz
    a = df.potential.a
    c = df.potential.c
    G = df.potential.G
    M = df.potential.M
    
    # todo
    
end