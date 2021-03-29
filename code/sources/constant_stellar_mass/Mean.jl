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
    halo = (a/b)^3*(1/q-1)*((1-xi)/2)^(-3/2)*(1+(r/b)^2)^(-3/2)
    disk = 1-eps/((1-x)^(2/3))
    return sqrt(1-x)*((1-xi)/2)^(3/4)*sqrt(bulb+halo+disk)
end

function dOmegadxi(xi::Float64,x::Float64=0.0, q::Float64=1.0)
    rc = a/c
    rb = a/b
    num = (-3+3*eps/(1-x)^(2/3)+(12*(-1+q)*rb^5*sqrt((2-2*xi)/(1-xi+rb^2*(1+xi))))/(q*sqrt(1-xi)*(1-xi+rb^2*(1+xi))^2)
           +(12*rc^5*x*sqrt((2-2*xi)/(1-xi+rc^2*(1+xi)))/((-1+x)*sqrt(1-xi)*(1-xi+rc^2*(1+xi))^2)))
    den = (4*2^(3/4)*(1-xi)^(1/4)*sqrt(1-eps/(1-x)^(2/3)+(2*sqrt(2)*(-1+1/q)*rb^3*(1-(rb^2*(1+xi))/(-1+xi))^(-3/2))
            /(1-xi)^(3/2)+(2*sqrt(2)*x*rc^3*(1-(rc^2*(1+xi))/(-1+xi))^(-3/2))
                    /((1-x)*(1-xi)^(3/2))))

    return sqrt(1-x)*(num/den)
end

function kappaSqOverTwoOmega(xi::Float64,x::Float64=0.0, q::Float64=1.0)
    omega = Omega(xi,x,q)
    domega = dOmegadxi(xi,x,q)
#    return (1-eps)^(1/2)*((1-xi)/2)^(3/4)*(1/2+(3/2)*(1-xi)/2)
    return 2*omega*(1+(1+xi)*(1-xi)*domega/(2*omega))
end

##############################
# Square Distance between two vectors
##############################

function dist2(v1::Array{Complex{Float64},1}, v2::Array{Complex{Float64},1})
    sizev1 = length(v1)
    sizev2 = length(v2)
    SqNorm = 0.0
    if (sizev1 <= sizev2)
        for k=1:sizev1
            SqNorm += abs2(v1[k]-v2[k])
        end
        # for k=sizev1+1:sizev2
        #     SqNorm += abs2(v2[k])
        # end
    else
        for k=1:sizev2
            SqNorm += abs2(v1[k]-v2[k])
        end
        # for k=sizev2+1:sizev1
        #     SqNorm += abs2(v1[k])
        # end
    end
    return SqNorm
end

##############################
# Test
##############################

function dOmegadxiNum(xi::Float64,x::Float64=0.0, q::Float64=1.0,prec::Float64=0.0001)
    return (Omega(xi+prec,x,q)-Omega(xi-prec,x,q))/(2*prec)
end
