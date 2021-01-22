#####################################
# Matrix elements
#####################################

using QuadGK
using Calculus

function _I(a::Float64, m::Int, l::Int, n::Int)
    if (l>=m && n>=m)
        local integrand
        let integrand
        integrand = (xi->((1-xi)/2)^a*redLegendrePol(xi,l,m)*redLegendrePol(xi,n,m))
        return quadgk(integrand,-1,1)[1]
    else
        return 0.0
    end
end

function _J(a::Float64, m::Int, l::Int, n::Int)
    if (l>=m && n>=m)
        local integrand
        let integrand
        integrand = (xi->((1-xi)/2)^a*((1+xi)/2)^(-1)*redLegendrePol(xi,l,m)*redLegendrePol(xi,n,m))
        return quadgk(integrand,-1,1)[1]
    else
        return 0.0
    end
end

function _A(l::Int, n::Int, eps::Float64, m::Int=2)
    I = _I(3/4,m,l,n)
    return m*sqrt(1-eps)*I
end

function _B(l::Int, n::Int, eps::Float64, m::Int=2)
    Jm = _J(3/4,m,l-1,n)
    J  = _J(3/4,m,l,n)
    Jp = _J(3/4,m,l+1,n)
    return ((1/2)*(sqrt((2*l+1)*(l+m+1)*(l-m+1)/(2*l+3))*Jp + J
           - sqrt((2*l+1)*(l+m)*(l-m)/(2*l-1))*Jm))
end

function _C(l::Int, n::Int, eps::Float64, m::Int=2)
    J = _J(3/4,m,l,n)
    return m*J
end

function _D(l::Int, n::Int, eps::Float64, m::Int=2)
    Im = _J(3/4,m,l,n-1)
    I  = _J(3/4,m,l,n)
    Ip = _J(3/4,m,l,n+1)
    return ((1/2)*(1/(2*n+1)-eps/3)*(sqrt((2*n+1)*(n+m+1)*(n-m+1)/(2*n+3))*Ip 
           - I + sqrt((2*n+1)*(n+m)*(n-m)/(2*n-1))*Im))
end

function _F(l::Int, n::Int, eps::Float64, m::Int=2)
    I = _I(3/4,m,l,n)
    return 2*sqrt(1-eps)*I
end

function _G(l::Int, n::Int, eps::Float64, m::Int=2)
    I = _I(3/4,m,l,n)
    return -m*(1/(2*n+1)-eps/3)*I
end

function _H(l::Int, n::Int, eps::Float64, m::Int=2)
    Im = _I(3/4,m,l,n)
    Ip = _I(7/4,m,l,n)
    return (1/2)*sqrt(1-eps)*(Im+3*Ip)
end
