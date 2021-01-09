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
# https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Generalization_via_hypergeometric_functions

function LegendrePol(z::Float64, n::Int64, m::Int64)
    # use a recursion identity
    # stock previous polynomial coefficients
    # take into argument nbBasis
    # maybe make nbBasis an command-line argument to be parsed?
    
    # recursion formula: (n-m+1)P_{n+1}^m(x) = (2n+1)xP_n^m(x)-(n+m)P_{n-1}^n(x)
    if (n < m)
        # P = 0
    else if (n == m)
        # P_n^n(x) = (-1)^n*(2n-1)!! * (1-x^2)^(n/2)
        # recursion on n previous n?
        # P_{n+1}^{n+1}(x) = -(2n+1)*sqrt(1-x^2)*P_n^n(x)
    else if (n == m+1)
        # P_{n+1}^n(x) = x*(2n+1)*P_n^n(x)
    else
        # (n-m+1)P_{n+1}^m(x) = (2n+1)xP_n^m(x)-(n+m)P_{n-1}^n(x)
        # hence
        # P_{n+1}^m(x) = (2n+1)/(n-m+1) *xP_n^m(x)-(n+m)/(n-m+1) *P_{n-1}^n(x)
end

#P_n^|m|(x)
function redLegendrePol(x::Float64, n::Int64, m::Int64)
    local norm
    let norm
    norm = sqrt(gamma(n-abs(m)+1)*(2*n+1)/(2*gamma(n+abs(m) ) ))
    return norm*LegendrePol(x,n,abs(m))
    end
end
