#####
#
#   Structures
#
#####

abstract type ShuDF <: DiscDF end

struct ExpShuDisk{modelT<:TaperedMestel} <: ShuDF
    potential::modelT       # potential model
    Rd::Float64             # Exponential disk scale factor
    G::Float64              # gravitational constant
    M::Float64              # Mass of the stellar disk
    Q::Float64              # Toomre parameter
    Lzmin::Float64          # Angular momentum inner taper
    Lzmax::Float64          # Angular momentum outer taper
    Lc::Float64             # Angular momentum smoothing scale
    isOdd::Bool
end

struct ExpShuDiskEven{modelT<:TaperedMestel} <: ShuDF
    potential::modelT       # potential model
    df_bare::ExpShuDisk     # Bare DF of the Shu exponential disk
    isOdd::Bool
end

struct ExpShuDiskOdd{modelT<:TaperedMestel} <: ShuDF
    potential::modelT       # potential model
    df_bare::ExpShuDisk     # Bare DF of the Shu exponential disk
    isOdd::Bool
end

#####
#
#   Constructors
#
#####

"""
ExpShuDiskEven([potential])

Even component of the Shu exponential disc distribution function. 
"""
function ExpShuDiskEven(;potential::TaperedMestel=TaperedMestel(0.5, 0.9, 1.0), Rd::Float64=1.0, G::Float64=1.0,
                    M::Float64=1.0, Q::Float64=1.2, Lzmin::Float64=0.0, Lzmax::Float64=5.381346952, Lc::Float64=0.1)
    return ExpShuDiskEven(potential, ExpShuDisk(potential, Rd, G, M, Q, Lzmin, Lzmax, Lc, false), false)
end

"""
ExpShuDiskEven([potential])

Odd component of the Shu exponential disc distribution function. 
"""
function ExpShuDiskOdd(;potential::TaperedMestel=TaperedMestel(0.5, 0.9, 1.0), Rd::Float64=1.0, G::Float64=1.0,
                    M::Float64=1.0, Q::Float64=1.2, Lzmin::Float64=0.0, Lzmax::Float64=5.381346952, Lc::Float64=0.1)
    return ExpShuDiskOdd(potential, ExpShuDisk(potential, Rd, G, M, Q, Lzmin, Lzmax, Lc, false), true)
end


#####
#
#   Auxiliary functions
#
#####

# Interpolation table
const tab_ck = [
    2.7515500106E+00,
    1.5649678643E+00,
    1.4106149953E+00,
    4.2634681533E-01,
    4.6164960163E-02,
    -5.9820363335E-01,
    -9.0880169334E-01,
    -1.2340267366E+00,
    -1.2574009691E+00,
    -1.2454080880E+00,
    -9.9301470728E-01,
    -7.7158623586E-01,
    -3.9922325904E-01,
    -1.4214586548E-01,
    1.9021152764E-01,
    3.5628862601E-01,
    5.5669823137E-01,
    5.8212341681E-01,
    6.3290433779E-01,
    5.3215838631E-01,
    4.7099677903E-01,
    3.0110453994E-01,
    1.9509340818E-01,
    2.5392008673E-02,
    -6.1548703184E-02,
    -1.7553072991E-01,
    -2.0610794749E-01,
    -2.4800897455E-01,
    -2.1856517044E-01,
    -2.0324543635E-01,
    -1.3765050455E-01,
    -9.7312779206E-02,
    -3.3361017144E-02,
    5.8501600222E-03,
    4.1967498659E-02,
    5.5072797087E-02,
    6.3937120145E-02,
    5.2623590236E-02,
    4.7939319156E-02,
    2.3340523655E-02,
]

function _bar_Lz(Lz::Float64, df::ExpShuDisk)
    Lzmin = df.Lzmin
    Lzmax = df.Lzmax
    return min(max(abs(Lz), Lzmin), Lzmax)
end

function _omega_shu(x::Float64)
    if (x <= -1.0)
        return 0.0
    elseif (x < 1.0)
        return 0.5 + 0.75*x - 0.25*x^3
    else
        return 1.0
    end
end

function _lambda(barLz::Float64, df::ExpShuDisk)
    Lzmax = df.Lzmax
    return (2*barLz-Lzmax)/Lzmax
end

function _g(barLz::Float64, df::ExpShuDisk)
    sum = 0.5*tab_ck_shu[1]
    lambda = _lambda(barLz, df)
    N = length(tab_ck_shu)

    for k=2:N 
        cheb = cos((k-1)*acos(lambda)) # T_{k-1}[lambda]
        sum += tab_ck_shu[k] * cheb
    end
    return sum 
end

function _Rc(Lz::Float64, potential::TaperedMestel)
    a = potential.R0
    V0 = potential.V0
    Rc2 = (Lz^2 + sqrt(Lz^4 + 4.0*a^2*V0^2*Lz^2))/(2.0*V0^2)
    return sqrt(Rc2)
end

function _Ec(Lz::Float64, potential::TaperedMestel)
    if (Lz > 0)
        return OrbitalElements._ψeff(_Rc(Lz, potential), Lz, potential)
    else
        return OrbitalElements.ψ(0.0, potential)
    end
end

function _Sigma(R::Float64, df::ExpShuDisk)
    M = df.M
    Rd = df.Rd
    # return M/(2*pi*Rd^2) * exp(-R/Rd)
    return 0.97 * M/(2*pi*Rd^2) * exp(-R/Rd) # Active fraction = 0.97
end

function _kappa(R::Float64, potential::TaperedMestel)
    a = potential.R0
    V0 = potential.V0
    return sqrt(2.0)*V0*sqrt(R^2+2*a^2)/(R^2+a^2)
end

function _sigma_R(R::Float64, df::ExpShuDisk)
    potential = df.potential
    Q = df.Q
    G = df.G
    return 3.36*Q*G*_Sigma(R, df)/_kappa(R, potential)
end


function _omega_shu_gradx(x::Float64)
    if (abs(x) < 1.0)
        return 0.75 * (1.0 - x^2)
    else
        return 0.0
    end
end

function _g_gradbarLz(barLz::Float64, df::ExpShuDisk)
    Lzmax = df.Lzmax
    N = length(tab_ck_shu)

    sum = 0.0
    lambda = _lambda(barLz, df)
    for k=2:N 
        if (-1.0 < lambda < 1.0)
            cheb_grad = (k-1)*sin((k-1)*acos(lambda))/sqrt(1.0-lambda^2) # T'_{k-1}[lambda]
        elseif (lambda == -1.0)
            cheb_grad = (-1)^k * k^2 # T'_{k-1}[-1]
        else
            cheb_grad = k^2 # T'_{k-1}[1]
        end
        sum += tab_ck_shu[k] * cheb_grad
    end
    return 2.0/Lzmax * sum 
end

function _bar_Lz_gradLz(Lz::Float64, df::ExpShuDisk)
    Lzmin = df.Lzmin
    Lzmax = df.Lzmax

    if (Lzmin < abs(Lz) < Lzmax)
        return sign(Lz)
    else
        return 0.0
    end
end

function _Rc_gradLz(Lz::Float64, potential::TaperedMestel)
    V0 = potential.V0
    a = potential.R0

    num = Lz^2 + 4*a^2*V0^2 + sqrt(Lz^4 + 4*a^2*V0^2*Lz^2)
    den = Lz^3 + 4*a^2*V0^2*Lz 
    Rc = _Rc(Lz, potential)
    return 0.5*Rc*num/den
end

function _Ec_gradLz(Lz::Float64, potential::TaperedMestel)
    V0 = potential.V0
    a = potential.R0

    if (Lz > 0)
        Rc = _Rc(Lz, potential)
        return sqrt(OrbitalElements.dψ(Rc, potential)/Rc)
    else
        return V0/a
    end
end

function _kappa_gradR(R::Float64, potential::TaperedMestel)
    V0 = potential.V0
    a = potential.R0

    return -(2*V0^2*R*(R^2+3*a^2))/(_kappa(R, potential)*(R^2+a^2)^3)
end

function _sigma_R_gradR(R::Float64, df::ExpShuDisk)
    Rd = df.Rd
    potential = df.potential

    return -_sigma_R(R,df) * (1.0/Rd + _kappa_gradR(R, potential)/_kappa(R, potential))
end


#####
#
#   Function and gradients for the bare DF
#
#####


"""
    ExpShuDisk(EL::Tuple{Float64,Float64},df::ExpShuDisk)
Shu exponential disc distribution function.
"""
function DistributionFunction(EL::Tuple{Float64,Float64},df::ExpShuDisk)::Float64
    E, Lz = EL
    Lc = df.Lc
    potential = df.potential

    barLz = _bar_Lz(Lz, df)
    omegaLz = _omega_shu(Lz/Lc)
    gLz = _g(barLz, df)
    Ec = _Ec(barLz, potential)
    Rc = _Rc(barLz, potential)
    sigma = _sigma_R(Rc, df)
    expELz = exp(-(E-Ec)/(sigma^2))

    if (E > Ec)
        return omegaLz * gLz * expELz
    else
        return 0.0
    end

end

"""
    DFDE(EL::Tuple{Float64,Float64},df::ExpShuDisk)
Shu exponential disc DF derivative w.r.t. E.
"""
function DFDE(EL::Tuple{Float64,Float64},df::ExpShuDisk)::Float64
    E, Lz = EL
    potential = df.potential

    barLz = _bar_Lz(Lz, df)
    Rc = _Rc(barLz, potential)
    sigma = _sigma_R(Rc, df)

    return -DistributionFunction(EL, df)/sigma^2
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::ExpShuDisk)
Shu exponential disc DF derivative w.r.t. L.
"""
function DFDL(EL::Tuple{Float64,Float64},df::ExpShuDisk)::Float64
    E, Lz = EL
    Lc = df.Lc
    potential = df.potential

    barLz = _bar_Lz(Lz, df)
    omegaLz = _omega_shu(Lz/Lc)
    gLz = _g(barLz, df)
    Ec = _Ec(barLz, potential)
    Rc = _Rc(barLz, potential)
    sigma = _sigma_R(Rc, df)
    expELz = exp(-(E-Ec)/(sigma^2))

    omegap = _omega_shu_gradx(Lz/Lc)
    gp = _g_gradbarLz(barLz, df)
    barLz_grad = _bar_Lz_gradLz(Lz, df)
    Ecp = _Ec_gradLz(barLz, potential)
    Rcp = _Rc_gradLz(barLz, potential)
    sigmap = _sigma_R_gradR(Rc, df)

    if (E <= Ec)
        return 0.0
    else

        bracket = Ecp/sigma^2 + 2*(E-Ec)*Rcp*sigmap/sigma^3
        DF = DistributionFunction(EL, df)

        term1 = omegap/Lc * gLz * expELz
        term2 = omegaLz * barLz_grad * gp * expELz
        term3 = barLz_grad * bracket * DF

        return term1 + term2 + term3
    end
end


#####
#
#   Function and gradients for the even DF
#
#####

"""
    ExpShuDisk(EL::Tuple{Float64,Float64},df::ExpShuDiskEven)
Shu exponential disc distribution function (even component).
"""
function DistributionFunction(EL::Tuple{Float64,Float64},df::ExpShuDiskEven)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DistributionFunction((E,L), df_bare) + DistributionFunction((E,-L), df_bare))
end

"""
    DFDE(EL::Tuple{Float64,Float64},df::ExpShuDiskEven)
Shu exponential disc DF derivative w.r.t. E (even component).
"""
function DFDE(EL::Tuple{Float64,Float64},df::ExpShuDiskEven)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DFDE((E,L), df_bare) + DFDE((E,-L), df_bare))
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::ExpShuDiskEven)
Shu exponential disc DF derivative w.r.t. L (even component).
"""
function DFDL(EL::Tuple{Float64,Float64},df::ExpShuDiskEven)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DFDL((E,L), df_bare) - DFDL((E,-L), df_bare))
end


#####
#
#   Function and gradients for the odd DF
#
#####

"""
    ExpShuDisk(EL::Tuple{Float64,Float64},df::ExpShuDiskOdd)
Shu exponential disc distribution function (odd component).
"""
function DistributionFunction(EL::Tuple{Float64,Float64},df::ExpShuDiskOdd)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DistributionFunction((E,L), df_bare) - DistributionFunction((E,-L), df_bare))
end

"""
    DFDE(EL::Tuple{Float64,Float64},df::ExpShuDiskOdd)
Shu exponential disc DF derivative w.r.t. E (odd component).
"""
function DFDE(EL::Tuple{Float64,Float64},df::ExpShuDiskOdd)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DFDE((E,L), df_bare) - DFDE((E,-L), df_bare))
end

"""
    dFdL(EL::Tuple{Float64,Float64},df::ExpShuDiskOdd)
Shu exponential disc DF derivative w.r.t. L (odd component).
"""
function DFDL(EL::Tuple{Float64,Float64},df::ExpShuDiskOdd)::Float64
    df_bare = df.df_bare
    E, L = EL
    return 0.5 * (DFDL((E,L), df_bare) + DFDL((E,-L), df_bare))
end