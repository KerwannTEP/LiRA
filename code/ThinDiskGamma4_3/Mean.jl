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
end

# Legendre polynomial through _₂F₁
# https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Generalization_via_hypergeometric_functions

function LegendrePol(x::Float64, n::Int, m::Int)
    # use a recursion identity
    # stock previous polynomial coefficients
    # take into argument nbBasis
    # maybe make nbBasis an command-line argument to be parsed?
 
    # recursion formula: (n-m+1)P_{n+1}^m(x) = (2n+1)xP_n^m(x)-(n+m)P_{n-1}^n(x)
    if (n < m)
        return 0.0
        # P = 0
    elseif (n == m)
        # P_n^n(x) = (-1)^n*(2n-1)!! * (1-x^2)^(n/2)
        # recursion on n previous n?
        # P_{n+1}^{n+1}(x) = -(2n+1)*sqrt(1-x^2)*P_n^n(x)
        if (n == 0)
            return 1.0
        else
            return -(2*n-1)*sqrt(1-x^2)*LegendrePol(x,n-1,n-1)
        end
    elseif (n == m+1)
        # P_{m+1}^m(x) = x*(2m+1)*P_m^m(x)
        return x*(2*m+1)*LegendrePol(x,m,m)
    else
        # (n-m+1)P_{n+1}^m(x) = (2n+1)xP_n^m(x)-(n+m)P_{n-1}^n(x)
        # hence
        # P_{n+1}^m(x) = (2n+1)/(n-m+1) *xP_n^m(x)-(n+m)/(n-m+1) *P_{n-1}^m(x)
        return (2*n-1)/(n-m) *x*LegendrePol(x,n-1,m) - (n-1+m)/(n-m) *LegendrePol(x,n-2,m)
    end
end

#P_n^|m|(x)
function redLegendrePol(x::Float64, n::Int, m::Int)
    local norm
    let norm
    norm = sqrt(gamma(n-abs(m)+1)*(2*n+1)/(2*gamma(n+abs(m) ) ))
    return norm*LegendrePol(x,n,abs(m))
    end
end
