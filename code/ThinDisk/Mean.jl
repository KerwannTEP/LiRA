#####################################
# Useful functions
#####################################

function _xi(r::Float64, a::Float64)
    return (r^2 - a^2)/(r^2 + a^2)
end

# r/a ; xi in ]-1,1[
function redRadius(xi::Float64)
    return sqrt((1+xi)/(1-xi))

function redRadialVelocity0(xi::Float64, eps::Float64)
    return 0.0
end

function redTangentVelocity0(xi::Float64, eps::Float64)
    return ((1-xi)/2)^(1/4)*((1+xi)/2)^(1/2)*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 - 2))
end

function redTotPotential(xi::Float64, eps::Float64)
    return -((1-xi)/2)^(1/2) + eps/(3*(polyIndex-1)) * ((1-xi)/2)^(3*(polyIndex-1)/2)
end

function redGravPotential(xi::Float64)
    return -((1-xi)/2)^(1/2)
end

function redDensity(xi::Float64)
    return ((1-xi)/2)^(3/2)
end

function redAngVelocity(xi::Float64, eps::Float64)
    return redTangentVelocity(xi,eps) / redRadius(xi)
end

function dRedAngVelocitydr(xi::Float64, eps::Float64)
    redRad = redRadius(xi)
    dxidr = 4*(redRad)/((redRad)^2+1)^2
    dvthdxi = (1/4)*(-1/2)*((1-xi)/2)^(-3/4)*((1+xi)/2)^(1/2)*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 - 2))
            + ((1-xi)/2)^(1/4)*(1/2)*(1/2)*((1+xi)/2)^(-1/2)*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 - 2))
            + ((1-xi)/2)^(1/4)*((1+xi)/2)^(1/2)*(-eps)*(3*polyIndex/2 - 2)*(-1/2)*(1-eps*((1-xi)/2)^(3*polyIndex/2 - 3)) / (2*sqrt(1-eps*((1-xi)/2)^(3*polyIndex/2 - 2)))
    return (dxidr*dvthdxi*redRad-redTangentVelocity0(xi,eps))/(redRad)^2
end

function redSqEpiFreqOverFourOmegaSq(xi::Float64, eps::Float64)
    omega = redAngVelocity(xi,eps)
    return (1+redRadius(xi)/(2*omega) * dxi/dr * dRedAngVelocitydxi(xi,eps))
end
