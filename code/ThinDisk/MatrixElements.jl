#####################################
# Matrix elements
#####################################

using QuadGK

function _A(l::Int64, n::Int64, eps::Float64, m::Int64=2)
    local integrand
    let integrand
    integrand = (xi->redLegendrePol(xi,l,m)*redAngularFreq(xi,eps)*redLegendrePol(xi,n,m))
    return abs(m)*quadgk(integrand,-1,1)[1]
    end
end

