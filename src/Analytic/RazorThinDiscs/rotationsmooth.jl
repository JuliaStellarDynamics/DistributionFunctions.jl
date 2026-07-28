
"""
ToomreDiscSmoothOdd([potential])

Toomre disc distribution function.
"""
function ToomreDiscSmoothOdd(;df_even::ToomreDisc=ToomreDisc(),a::Float64=1.0,Lc::Float64=1.0)
    @assert (0.0 < a < 1.0) "ERROR : a should be strictly between 0.0 and 1.0"
    potential = df_even.potential
    G = df_even.G
    return ToomreDiscSmoothOdd(potential,df_even,a,Lc,G,true)
end


function _gRot(x::Float64)
    if (x <= -1.0)
        return -1.0
    elseif (x < 1.0)
        return x * (1.5 - 0.5*x*x)
    else
        return 1.0
    end
end

function _gRot_a(x::Float64, a::Float64)
    return _gRot(x/a)
end

function _d_gRot_dx(x::Float64)
    if abs(x) >= 1.0
        return 0.0
    else
        return 1.5 * (1.0 - x*x)
    end
end

function _d_gRot_a_dx(x::Float64, a::Float64)
    return _d_gRot_dx(x/a)/a
end

"""
    DistributionFunction(EL::Tuple{Float64,Float64},df::ToomreDiscSmoothOdd)
Toomre distribution function with smooth rotation (odd component).
"""
function DistributionFunction(EL::Tuple{Float64,Float64},df::ToomreDiscSmoothOdd)::Float64
    df_even = df.df_even
    Lc = df.Lc
    a = df.a
    E,L = EL
    x = L/Lc
    return DistributionFunction(EL,df_even) * _gRot_a(x, a)
end

"""
    dFdE(EL::Tuple{Float64,Float64},df::ToomreDiscSmoothOdd)
Toomre DF derivative w.r.t. E (smooth odd component).
"""
function DFDE(EL::Tuple{Float64,Float64},df::ToomreDiscSmoothOdd)::Float64
    df_even = df.df_even
    Lc = df.Lc
    a = df.a
    E,L = EL
    x = L/Lc
    return DFDE(EL,df_even) * _gRot_a(x, a)
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::ToomreDiscSmoothOdd)
Toomre DF derivative w.r.t. L (smooth odd component).
"""
function DFDL(EL::Tuple{Float64,Float64},df::ToomreDiscSmoothOdd)::Float64
    df_even = df.df_even
    Lc = df.Lc
    a = df.a
    E,L = EL
    x = L/Lc
    return DFDL(EL,df_even) * _gRot_a(x, a) + DistributionFunction(EL,df_even) * _d_gRot_a_dx(x, a)/Lc
end