function _xi(r::Float64)
    r_a = r/a
    return (r_a^2-1)/(r_a^2+1)
end

function radius(xi::Float64)
    return a*sqrt((1+xi)/(1-xi))
end

function _x(frac_bulb::Float64=0.0)
    return frac_disk
end

function _q(frac_bulb::Float64=0.0, frac_DH::Float64=0.0)
    return (1-frac_bulb)/(frac_DH+1-frac_bulb)
end

function Omega(xi::Float64,x::Float64=0.0, q::Float64=1.0)
    r = radius(xi)
    bulb = (a/c)^3*x/(1-x)*((1-xi)/2)^(-3/2)*(1+(r/c)^2)^(-3/2)
    disk = 1/q-eps/((1-x)^(2/3))
#    return sqrt(1-eps)*((1-xi)/2)^(3/4)
    return sqrt(1-x)*((1-xi)/2)^(3/4)*sqrt(bulb+disk)
end

function dOmegadxi(xi::Float64,x::Float64=0.0, q::Float64=1.0)
    r = a/c
    xq = x*q

    num = 3*((-1+x+q*(1-x)^(1/3)*eps)*(1-xi)^(5/2)
            +r^4*(-1+x+q*(1-x)^(1/3)*eps)*sqrt(1-xi)*(1+xi)^2
            -2*r^2*(-1+x+q*(1-x)^(1/3)*eps)*sqrt(1-xi)*(-1+xi^2)
            -4*q*r^5*x*sqrt((2-2*xi)/(1-xi+r^2*(1+xi))))

    den = (4*q*(-1+x)*(2-2*xi)^(3/4)*(1-xi+r^2*(1+xi))^2*sqrt(1/q
          - eps/(1-x)^(2/3) - r^3*x*((2-2*xi)/(1-xi+r^2*(1+xi)))^(3/2)
          /((-1+x)*(1-xi)^(3/2))))

    return sqrt(1-x)*(-num/den)
end

#correct this
function kappaSqOverTwoOmega(xi::Float64,x::Float64=0.0, q::Float64=1.0)
    omega = Omega(xi,x,q)
    domega = dOmegadxi(xi,x,q)
#    return (1-eps)^(1/2)*((1-xi)/2)^(3/4)*(1/2+(3/2)*(1-xi)/2)
    return 2*omega*(1+(1+xi)*(1-xi)*domega/(2*omega))
end

##############################
# Test
##############################

function dOmegadxiNum(xi::Float64,x::Float64=0.0, q::Float64=1.0,prec::Float64=0.0001)
    return (Omega(xi+prec,x,q)-Omega(xi-prec,x,q))/(2*prec)
end
