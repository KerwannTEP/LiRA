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
        tabEigValsMln[k] = sqrt(x/q)*eigvals(tabTruncMln[k])
    end
end

function getPhysicalEigvals()
    liEgv = zeros(Complex{Float64},3*(N+1))
    liEigvals = tabEigValsMln[N+1]
    for iegv=1:3*(N+1) # take one of the eigenvalues of Mln(N+1)
        egv = liEigvals[iegv]
        liEgvCvCurrent = zeros(Complex{Float64},N)
        for k=1:N # compare to those of its truncatures
            dist_min = -1
            egvSelect = 0
            for kegv=1:3*k # at a given truncature, take the closest one
                egvTrunc = tabEigValsMln[k][kegv]
                dist = abs(egv,egvTrunc)
                if ((d<0) || (dist<dist_min))
                    egvSelect = egvTrunc
                    dist_min = dist
                end
            end
            if ((d>=0) && (dist_min<error_dis)) # check if the closest one is within the accepted error margin
                liEgvCvCurrent[k] = egvSelect
            else # if not, put Inf to discard
                liEgvCvCurrent[k] = Inf
            end

        nb_egv = 1 #(1 comes from initial separate egv at N+1)
        for k=1:N # count number of accepted truncature eigenvalues
            if (liEgvCvCurrent[k] != Inf)
                nb_egv += 1
            end
        end
        if (nb_egv >= threshold_nb_egv) # if enough eigenvalues 
            moy_egv = 0
            var_egv = 0

            for k=1:N # compute mean
                if (liEgvCvCurrent[k] != Inf)
                    moy_egv += liEgvCvCurrent[k]
                end
            end
            moy_egv = (moy_egv+egv)/nb_egv

            for k=1:N # compute variance
                if (liEgvCvCurrent[k] != Inf)
                    var_egv += (liEgvCvCurrent[k]-moy_egv)^2
                end
            end
            var_egv = (var_egv+(egv-moy_egv)^2)/nb_egv

            if (var_egv < error_var) # recover mean physical eigenvalue
                liEgv[iegv] = moy_egv
            else
                liEgv[iegv] = Inf
            end
        end
    end

    nb_phys_egv = 0
    for iegv=1:3*(N+1)
        if (liEgv[iegv] != Inf)
            nb_phys_egv += 1
        end
    end
    
    liPhysicalEgv = zeros(Complex{Float64},nb_phys_egv)
    index = 1
    for iegv=1:3*(N+1)
        if (liEgv[iegv] != Inf)
            liPhysicalEgv[index] = liEgv[iegv]
            index += 1
        end
    end

    return liPhysicalEgv
end


# Gets the element of the list li with maximum imaginary part
function get_max_Im(li) 
    len = length(li)
    max = -Inf
    egv = -Inf
    for k=1:len
        if (imag(li[k])>max)
            egv = li[k]
            max = imag(egv)
        end
    end
    return egv
end
