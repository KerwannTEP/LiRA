#####################################
# Useful functions
#####################################

using HypergeometricFunctions
using SpecialFunctions

function _xi(r::Float64, a::Float64)
    return (r^2 - a^2)/(r^2 + a^2)
end

# r/a ; xi in ]-1,1[
function redRadius(xi::Float64)
    return sqrt((1+xi)/(1-xi))

function redAngularFreq(xi::Float64, eps::Float64)
    return ((1-xi)/2)^(1/2)*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 -2))
end

function redDensity(xi::Float64)
    return ((1-xi)/2)^(3/2)
end

function redSqEpiFreqOverRedAngFreq(xi::Float64, eps::Float64)
    local fct, pref
    let fct, pref
    fct = derivative(x->(1-x)^(1/2)*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 -2)))
    pref = (1+xi)*(1-xi)^(1/2)/(4*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 -2)))
    return 2*((1-xi)/2)^(1/2)*(1+pref*fct(xi))
    end
end

# Legendre polynomial through _₂F₁
# https://functions.wolfram.com/Polynomials/LegendreP2/26/01/02/0004/

function LegendrePol(x::Float64, n::Int64, m::Int64)
    local pref
    let pref
    pref = gamma(2*n+1)/(gamma(n+1)*gamma(n-m+1))
    if (mod(n,2) == 0)
        pref *= 1/2^n
    else
        pref *= -1/2^n
    end
    return pref*(1-z)^(n-m/2)*(1+z)^(m/2)*_₂F₁(-n,m-n,-2*n,2/(1-z))
    end
end

#P_n^|m|(x)
function redLegendrePol(x::Float64, n::Int64, m::Int64)
    local norm
    let norm
    norm = sqrt(gamma(n-abs(m)+1)*(2*n+1)/(2*gamma(n+abs(m) ) ))
    return norm*LegendrePol(x,n,abs(m))
    end
end
