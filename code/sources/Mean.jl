function _xi(r::Float64)
    r_a = r/a
    return (r_a^2-1)/(r_a^2+1)
end

function radius(xi::Float64)
    return a*sqrt((1+xi)/(1-xi))
end

function Omega(xi::Float64,x::Float64=1.0, q::Float64=1.0)
    r = radius(xi)
    bulb = (a/c)^3*q*(1-x)/x*((1-xi)/2)^(-3/2)*(1+(r/c)^2)^(-3/2)
    disk = 1-eps*(q/x)^(2/3)
    return ((1-xi)/2)^(3/4)*sqrt(bulb+disk)
end

function dOmegadxi(xi::Float64,x::Float64=1.0, q::Float64=1.0)
    r = a/c
    xq = x/q
    num = 3*(-4*q*r^5*sqrt((2-2*xi)/(1-xi+r^2*(1+xi)))
        + x*((-1+xq^(2/3)*eps)*(1-xi)^(5/2)
        + r^4*(-1+xq^(2/3)*eps)*sqrt(1-xi)*(1+xi)^2
        - 2*r^2*(-1+xq^(2/3)*eps)*sqrt(1-xi)*(-1+xi^2)
        + 4*q*r^5*sqrt((2-2*xi)/(1-xi+r^2*(1+xi)))))
    den = (4*x*(2-2*xi)^(3/4)*(1-xi+r^2*(1+xi))^2
          *sqrt(1-xq^(2/3)*eps+(2*sqrt(2)*q*r^3*(1-x)
          /(1+r^2*(1+xi)/(1-xi))^(3/2))
          /(x*(1-xi)^(3/2))))
    return num/den
end

function alphaSqOverTwoOmega(xi::Float64,x::Float64=1.0, q::Float64=1.0)
    omega = Omega(xi,x,q)
    domega = dOmegadxi(xi,x,q)
    return 2*omega*(1+(1+xi)*(1-xi)*domega/(2*omega))
end

##############################
# Test
##############################

function dOmegadxiNum(xi::Float64,prec::Float64=0.0001,x::Float64=1.0, q::Float64=1.0)
    return (Omega(xi+eps,x,q)-Omega(xi-eps,x,q))/(2*eps)
end
