#####################################
# Useful functions
#####################################

using Calculus
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

function RedSqEpiFreqOverRedAngFreq(xi::Float64, eps::Float64)
    local fct, pref
    let fct, pref
    fct = derivative(x->(1-x)^(1/2)*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 -2)))
    pref = (1+xi)*(1-xi)^(1/2)/(4*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 -2)))
    return 2*((1-xi)/2)^(1/2)*(1+pref*fct(xi))
    end
end

# Legendre polynomial through _₂F₁
# P_n^m(z)= 1/gamma(1-m) * ((1+z)/(1-z))^(m/2) *_₂F₁(-n,n+1,1-m,(1-z)/2) (finite sum)

# P_n^|m|
function redLegendrePol(xi::Float64, n::Int64, m::Int64)
    local norm
    let norm
    norm = sqrt( gamma(n-abs(m)+1)*(2*n+1)/(2*gamma(n+abs(m)+1)) )
    return norm * 1/gamma(1-abs(m)) * ((1+xi)/(1-xi))^(m/2) *_₂F₁(-n,n+1,1-m,(1-xi)/2)
    end
end
