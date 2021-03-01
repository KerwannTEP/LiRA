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

function alphaSqOverTwoOmega(xi::Float64)
    return 1 # write expression
end
