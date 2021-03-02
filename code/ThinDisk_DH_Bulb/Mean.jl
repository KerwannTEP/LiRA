function _xi(r::Float64)
    r_a = r/a
    return (r_a^2-1)/(r_a^2+1)
end

function radius(xi::Float64)
    return a*sqrt((1+xi)/(1-xi))
end

function Omega(xi::Float64)
    bulb = (a/c)^3*q*(1-x)/x*((1-xi)/2)^(-3/2)*(1+(r/c)^2)^(-3/2)
    disk = 1-eps*(q/x)^(2/3)
    return ((1-xi)/2)^(3/4)*sqrt(bulb+disk)
end

function dOmegadxi(xi::Float64)
    r = a/c
    xq = x/q
    num = (-4*q*r^5*sqrt((2-2*xi)/(1-xi+r^2*(1+xi)))
        + x*((-1+xq^(2/3)*eps)*(1-xi)^(5/2)
        + r^4*(-1+xq^(2/3)*eps)*sqrt(1-xi)*(1+xi)^2
        - 2*r^2*(-1+xq^(2/3)*eps)*sqrt(1-xi)*(1-xi^2)
        + 4*q*r^5*sqrt((2-2*xi)/(1-xi+r^2*(1+xi)))))
    den = (4*x*(2-2*xi)^(3/4)*(1-xi+r^2*(1+xi))^2
          *sqrt(1-xq^(2/3)*eps+(2*sqrt(2)*q*r^3*(1-x)
          /(1+r^2*(1+xi)/(1-xi))^(3/2)))
          /(x*(1-xi)^(3/2))))
    return num/den
end

function alphaSqOverTwoOmega(xi::Float64)
    omega = Omega(xi)
    domega = dOmegadxi(xi)
    return 2*omega*(1+(1+xi)*(1-xi)*domega/(2*omega)) 
end
