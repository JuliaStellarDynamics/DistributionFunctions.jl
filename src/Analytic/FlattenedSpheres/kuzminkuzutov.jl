

"""
KuzminKuzutov([potential])

Isotropic Plummer distribution function. Uses OrbitalElements.NumericalPlummer by default.
"""
function KuzminKuzutovELz(;potential::KuzminKuzutovPotential=AnalyticKuzminKuzutov())
    return KuzminKuzutovELz(potential)
end



"""
The distribution function for the isotropic Plummer model

Using the integral definition from Batsleer & Dejonghe (1993)
https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
"""
function DistributionFunction(ELz::Tuple{Float64,Float64}, df::KuzminKuzutovELz)
    E,Lz = ELz

    scaleEnergy = _tE(E)
    scaleLz = _tLz(Lz)

    # tc = c/(a+c)
    # ta = a/(a+c)
    scaleDF = df.potential.M/(df.potential.G*df.potential.M*(df.potential.a+df.potential.c))^(3/2)*(tc^2)/(2^(3/2)*pi^3*ta)
    z = 2*h*tE*tLz^2

    sum = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        num = (1-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        sum += num/den 
    end 

    sum *= f0*tE^(5/2)*1.0/nbK

    return sum 
end 

end


function x_eps(tE::Float64, eps::Float64, z::Float64, t::Float64)

    return 2.0*ta*tE*t*sqrt(1.0-t^2)/(1.0+eps*sqrt(z)*t)
end


# 
function F_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=nbK_default)

    tE = _tE(E)
    tLz = _tLz(Lz)
    f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2^(3/2)*pi^3*ta)
    z = 2*h*tE*tLz^2

    sum = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        num = (1-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        sum += num/den 
    end 

    sum *= f0*tE^(5/2)*1.0/nbK

    return sum 
end 

function F(E::Float64, Lz::Float64, nbK::Int64=100)

    return F_eps(-1.0,E,Lz,nbK) + F_eps(1.0,E,Lz,nbK)

end



function F_actions(Ju::Float64, Jv::Float64, Lz::Float64, nbK::Int64=nbK_default)

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    return F(E,Lz,nbK)

end

using Plots

# Derivatives

# For now, use the integral definition
# Batsleer & Dejonghe (1993)
# https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
function _dFdE_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=nbK_default)


    # dF/dtE
    tE = _tE(E)
    tLz = _tLz(Lz)
    f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2.0^(3/2)*pi^3*ta)
    z = 2.0*h*tE*tLz^2

    # Derivative on tE^2.5
    sum1 = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        num = (1.0-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1.0-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        sum1 += num/den 
    end 

    sum1 *= 2.5*f0*tE^(3/2)*1.0/nbK

    dzdE = 2.0*h*tLz^2


    sum2 = 0.0
    for k=1:nbK

        t = 1.0/nbK*(k-0.5)

        z = 2*h*tE*tLz^2
        xe = 2*ta*tE*t*sqrt(1-t^2)/(1+eps*sqrt(z)*t)

        dzdE = 2*h*tLz^2
        dxedE = 2*ta*t*sqrt(1-t^2)/(1+eps*sqrt(z)*t) - 2*ta*tE*t*sqrt(1-t^2)*(eps*0.5*dzdE/sqrt(z)*t)/(1+eps*sqrt(z)*t)^2


        num = (3+4*xe-xe^2)*(1-xe)*(1-t^2)+ 12*t^2
        den = (1-2*ta*tE*t*sqrt(1-t^2) + eps*sqrt(z)*t)^5

        dnumdE = (4*dxedE-2*dxedE*xe)*(1-xe)*(1-t^2)+  (3+4*xe-xe^2)*(-dxedE)*(1-t^2) 
        ddendE = 5*(1-2*ta*tE*t*sqrt(1-t^2) + eps*sqrt(z)*t)^4 * (-2*ta*t*sqrt(1-t^2) + eps*0.5*dzdE/sqrt(z)*t)

        sum2 += (1-t^2)*(dnumdE*den-num*ddendE)/den^2

    end 

    sum2 *= f0*tE^(5/2)*1.0/nbK

    sum = sum1 + sum2
    sum *= -(a+c)/(G*M)


    return sum 
end 




function _dFdE(E::Float64, Lz::Float64, nbK::Int64=nbK_default)

    return _dFdE_eps(-1.0,E,Lz,nbK) + _dFdE_eps(1.0,E,Lz,nbK)

end



# For now, use the integral definition
# Batsleer & Dejonghe (1993)
# https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
function _dFdLz_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=nbK_default)


    # dF/dtLz
    tE = _tE(E)
    tLz = _tLz(Lz)
    f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2.0^(3/2)*pi^3*ta)
    z = 2.0*h*tE*tLz^2

   
    dzdtLz = 2*h*tE*2*tLz

    # Derivative within the integral 
    sum2 = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        dxedtLz =  2*ta*tE*t*sqrt(1-t^2)*(-eps*t*0.5*dzdtLz/sqrt(z))/(1+eps*sqrt(z)*t)^2

        num = (1.0-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1.0-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        dnum = (1.0-t^2)*((4*dxedtLz-2*xe*dxedtLz)*(1-xe)*(1.0-t^2) + (3.0+4.0*xe-xe^2)*(-dxedtLz)*(1.0-t^2))
        dden = 5.0 * ( eps*t*(0.5*dzdtLz/sqrt(z)) ) * (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^4

        sum2 += (dnum*den-num*dden)/den^2 
    end 

    sum2 *= f0*tE^(5/2)*1.0/nbK


    sum2 *= 1/sqrt((a+c)*G*M)

    return sum2
end 

function _dFdLz(E::Float64, Lz::Float64, nbK::Int64=nbK_default)

    return _dFdLz_eps(-1.0,E,Lz,nbK) + _dFdLz_eps(1.0,E,Lz,nbK)

end

