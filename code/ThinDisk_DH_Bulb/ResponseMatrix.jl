using LinearAlgebra

const tabTruncMln = Vector{Array{Float64,2}}([zeros(Float64,3*k,3*k) for k=1:N+1])
const tabEigValsMln = Vector{Vector{Complex{Float64}}}([zeros(Float64,k) for k=1:N+1])

# Apply table_fill!() and matrix_fill!() beforehand.

function tabTruncMln!()
    for k=1:N+1
        for l=1:k
            for n=1:k
                tabTruncMln[k][l,n] = tabAln[l,n]
                tabTruncMln[k][l,k+n] = tabBln[l,n]
                tabTruncMln[k][l,2*k+n] = tabCln[l,n]
                tabTruncMln[k][k+l,n] = tabDln[l,n]
                tabTruncMln[k][k+l,k+n] = tabAln[l,n]
                tabTruncMln[k][k+l,2*k+n] = tabFln[l,n]
                tabTruncMln[k][2*k+l,n] = tabGln[l,n]
                tabTruncMln[k][2*k+l,k+n] = tabHln[l,n]
                tabTruncMln[k][2*k+l,2*k+n] = tabAln[l,n]
            end
        end
    end
end

function tabEigValsMln!()
    for k=1:N+1
        tabEigValsMln[k] = eigvals(tabTruncMln[k])
    end
end
