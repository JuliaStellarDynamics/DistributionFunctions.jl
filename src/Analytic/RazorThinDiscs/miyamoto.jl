
"""
ToomreDisc([potential])

Toomre disc distribution function.
"""
function ToomreDisc(;potential::ToomrePotential=NumericalToomre(),mM::Int64=1,G::Float64=1.0)
    return ToomreDisc(potential,mM,G,false)
end

"""
    ToomreDistribution(EL::Tuple{Float64,Float64},df::ToomreDisc)
Toomre distribution function.
"""
function DistributionFunction(EL::Tuple{Float64,Float64},df::ToomreDisc)::Float64

    return MiyamotoDistribution(EL,df)
end

"""
    dFdE(EL::Tuple{Float64,Float64},df::ToomreDisc)
Toomre DF derivative w.r.t. E.
"""
function DFDE(EL::Tuple{Float64,Float64},df::ToomreDisc)::Float64
    return MiyamotodFdE(EL,df)
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::ToomreDisc)
Toomre DF derivative w.r.t. E.
"""
function DFDL(EL::Tuple{Float64,Float64},df::ToomreDisc)::Float64
    return MiyamotodFdL(EL,df)
end





"""
ToomreDisc([potential])

Toomre disc distribution function.
"""
function ToomreDiscOdd(;potential::ToomrePotential=NumericalToomre(),mM::Int64=1,G::Float64=1.0)
    return ToomreDiscOdd(potential,mM,G,true)
end

"""
    ToomreDistribution(EL::Tuple{Float64,Float64},df::ToomreDisc)
Toomre distribution function.
"""
function DistributionFunction(EL::Tuple{Float64,Float64},df::ToomreDiscOdd)::Float64

    return MiyamotoDistributionOdd(EL,df)
end

"""
    dFdE(EL::Tuple{Float64,Float64},df::ToomreDisc)
Toomre DF derivative w.r.t. E.
"""
function DFDE(EL::Tuple{Float64,Float64},df::ToomreDiscOdd)::Float64
    return MiyamotodFdEOdd(EL,df)
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::ToomreDisc)
Toomre DF derivative w.r.t. E.
"""
function DFDL(EL::Tuple{Float64,Float64},df::ToomreDiscOdd)::Float64
    return MiyamotodFdLOdd(EL,df)
end