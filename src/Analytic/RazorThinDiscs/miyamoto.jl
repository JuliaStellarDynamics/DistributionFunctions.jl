using HypergeometricFunctions

#####
#
#   The Miyamoto DFs for Kuzmin-Toomre disc 
#
#####
"""
    miyamoto_distribution(E, L[, mM])

Miyamoto distribution function for Kuzmin-Toomre disc.
"""
function DistributionFunction(E::Float64, L::Float64; mM::Int64=1)
    return (
        (2mM + 3)
        * (-E)^(2mM + 2)
        * _₂F₁(-mM, -2-2mM, 1/2, -L^2/(2E))
        / (2 * pi^2)
    )
end

"""
    miyamoto_dFdE(E, L[, mM])

Miyamoto distribution derivative w.r.t. energy.
"""
function DFDE(E::Float64, L::Float64; mM::Int64=1)
    return (
        E^(2mM)
        * (1 + mM)
        * (3 + 2mM)
        * ( 
            L^2 * mM * _₂F₁(-1-2mM, 1-mM, 3/2, -L^2/(2E))
            + E * _₂F₁(-mM, -2-2mM, 1/2, -L^2/(2E))
        )
        / pi^2
    )
end

"""
    miyamoto_dFdL(E, L[, mM])

Miyamoto distribution derivative w.r.t. angular momentum.
"""
function DFDL(E::Float64, L::Float64; mM::Int64=1)

    return -(
        2
        * E^(1 + 2mM)
        * L
        * mM
        * (1 + mM)
        * (3 + 2mM)
        * _₂F₁(-1-2mM, 1-mM, 3/2, -L^2/(2E))
        / pi^2
    )
end

