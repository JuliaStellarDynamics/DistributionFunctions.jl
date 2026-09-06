struct TruncatedZangDiscEven{modelT<:MestelPotentials} <: ZangDF
    potential::modelT       # potential model
    df_bare::TruncatedZangDisc     # Bare DF of the Shu exponential disk
    isOdd::Bool
end

struct TruncatedZangDiscOdd{modelT<:MestelPotentials} <: ZangDF
    potential::modelT       # potential model
    df_bare::TruncatedZangDisc     # Bare DF of the Shu exponential disk
    isOdd::Bool
end

#####
#
#   Constructors
#
#####

function TruncatedZangDiscEven(;potential::MestelPotentials=MestelPotential(),q::IntorFloat=11.44,ν::Int64=4,Rin::Float64=1.0,μ::Int64=5,Rout::Float64=11.5,G::Float64=1.0)
    return ExpShuDiskEven(potential, ZangDisc(potential,q,ν,Rin,μ,Rout,G,false), false)
end

function TruncatedZangDiscOdd(;potential::MestelPotentials=MestelPotential(),q::IntorFloat=11.44,ν::Int64=4,Rin::Float64=1.0,μ::Int64=5,Rout::Float64=11.5,G::Float64=1.0)
    return ExpShuDiskOdd(potential, ZangDisc(potential,q,ν,Rin,μ,Rout,G,false), true)
end


#####
#
#   Function and gradients for the even DF
#
#####


function DistributionFunction(EL::Tuple{Float64,Float64},df::TruncatedZangDiscEven)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DistributionFunction((E,L), df_bare) + DistributionFunction((E,-L), df_bare))
end

function DFDE(EL::Tuple{Float64,Float64},df::TruncatedZangDiscEven)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DFDE((E,L), df_bare) + DFDE((E,-L), df_bare))
end

function DFDL(EL::Tuple{Float64,Float64},df::TruncatedZangDiscEven)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DFDL((E,L), df_bare) - DFDL((E,-L), df_bare))
end


#####
#
#   Function and gradients for the odd DF
#
#####


function DistributionFunction(EL::Tuple{Float64,Float64},df::TruncatedZangDiscOdd)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DistributionFunction((E,L), df_bare) - DistributionFunction((E,-L), df_bare))
end

function DFDE(EL::Tuple{Float64,Float64},df::TruncatedZangDiscOdd)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DFDE((E,L), df_bare) - DFDE((E,-L), df_bare))
end

function DFDL(EL::Tuple{Float64,Float64},df::TruncatedZangDiscOdd)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DFDL((E,L), df_bare) + DFDL((E,-L), df_bare))
end