#####################################
# Response matrix
#####################################

using LinearAlgebra

#       A B C
# M = ( D A F )
#       G H A
# iBasis from i=0 to nbBasis-1
function computeResponseMatrix(nbBasis::Int64, eps::Float64, m::Int64=2)
    M = zeros(3*nbBasis,3*nbBasis)
    elem = 0.0

    # make this computation parallel?

    for i=1:nbBasis
        for j=1:nbBasis

            # A component
            elem = _A(i-1,j-1,eps,m)
            M[i,j] = elem
            M[nbBasis+i,nbBasis+j] = elem
            M[2*nbBasis+i,2*nbBasis+j] = elem
        
            # B component
            elem = _B(i-1,j-1,eps,m)
            M[i,nbBasis+j] = elem

            # C component
            elem = _C(i-1,j-1,eps,m)
            M[i,2*nbBasis+j] = elem

            # D component
            elem = _D(i-1,j-1,eps,m)
            M[nbBasis+i,j] = elem

            # F component
            elem = _F(i-1,j-1,eps,m)
            M[nbBasis+i,2*nbBasis+j] = elem

            # G component
            elem = _G(i-1,j-1,eps,m)
            M[2*nbBasis+i,j] = elem

            # H component
            elem = _G(i-1,j-1,eps,m)
            M[2*nbBasis+i,nbBasis+j] = elem
        end
    end
    return M
end

function computeMaxEigenvalue(matrix::Array{Float64,2})
    D = eigvals(matrix)
    max = 0.0
    len = length(D)
    for i=1:len
        if abs(D[i])>max
            max = abs(D[i])
        end
    end
    return max
end

function growthRate(nbBasis::Int64, eps::Float64, m::Int64=2)
    responseMatrix = computeResponseMatrix(nbBasis,eps,m)
    growthRate = computeMaxEigenvalue(responseMatrix)
    return growthRate
end

# Select the true converging eigenvalues
