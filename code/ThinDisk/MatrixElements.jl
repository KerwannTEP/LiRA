#####################################
# Matrix elements
#####################################

using QuadGK
using Calculus

function _A(l::Int, n::Int, eps::Float64, m::Int=2)
    if (l>=m && n>=m)
        local integrand
        let integrand
        integrand = (xi->redLegendrePol(xi,l,m)*redAngularFreq(xi,eps)*redLegendrePol(xi,n,m))
        return abs(m)*quadgk(integrand,-1,1)[1]
        end
    else
        return 0.0
    end
end

function _B(l::Int, n::Int, eps::Float64, m::Int=2)
    if (l>=m && n>=m)
        local integrand, deriv
        let integrand, deriv
        deriv = derivative(x->((1-x)/2)^(-1/4) * redDensity(x)*redLegendrePol(x,n,m))
        integrand = (xi->redLegendrePol(xi,l,m)*((1-xi)/2)^(1/2)*deriv(xi))
        return 4*quadgk(integrand,-1,1)[1]
        end
    else
        return 0.0
    end
end

function _C(l::Int, n::Int, eps::Float64, m::Int=2)
    if (l>=m && n>=m)
        local integrand
        let integrand
        integrand = (xi->redLegendrePol(xi,l,m)*redDensity(xi)*((1+xi)/2)^(-1)
                     *((1-xi)/2)^(3/4)*redLegendrePol(xi,n,m))
        return abs(m)*quadgk(integrand,-1,1)[1]
        end
    else
        return 0.0
    end
end

function _D(l::Int, n::Int, eps::Float64, m::Int=2)
    if (l>=m && n>=m)
        local integrand, deriv
        let integrand, deriv
        deriv = derivative(x->(1/(2*n+1)-eps/3 *((1-x)/2)*(redDensity(x))^(polyIndex-2))
                *((1-x)/2)^(1/2)*redLegendrePol(x,n,m))
        integrand = (xi->redLegendrePol(xi,l,m)*((1-xi)/2)^(5/4)*((1+xi)/2)*deriv(xi))
        return 4*quadgk(integrand,-1,1)[1]
        end
    else
        return 0.0
    end
end

function _F(l::Int, n::Int, eps::Float64, m::Int=2)
    if (l>=m && n>=m)
        local integrand
        let integrand
        integrand = (xi->redLegendrePol(xi,l,m)*redAngularFreq(xi,eps)*redLegendrePol(xi,n,m))
        return 2*quadgk(integrand,-1,1)[1]
        end
    else
        return 0.0
    end
end

function _G(l::Int, n::Int, eps::Float64, m::Int=2)
    if (l>=m && n>=m)
        local integrand
        let integrand
        integrand = (xi->redLegendrePol(xi,l,m)*((1-xi)/2)^(3/4)*(1/(2*n+1)-eps/3*((1-xi)/2) 
                    * (redDensity(xi))^(polyIndex-2)) *redLegendrePol(xi,n,m)) 
        return -abs(m)*quadgk(integrand,-1,1)[1]
        end
    else
        return 0.0
    end
end

function _H(l::Int, n::Int, eps::Float64, m::Int=2)
    if (l>=m && n>=m)
        local integrand
        let integrand
        integrand = (xi->redLegendrePol(xi,l,m)*redSqEpiFreqOverRedAngFreq(xi,eps)
                     *redLegendrePol(xi,n,m)) 
        return quadgk(integrand,-1,1)[1]
        end
    else
        return 0.0
    end
end
