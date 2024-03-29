"""
    tabTruncMln_serial

Table containing the truncated response matrices from size `3*1` to `3*N`.
"""
tabTruncMln_serial = Vector{Array{Float64,2}}([zeros(Float64,3*k,3*k) for k=1:N+1])

"""
    tabEigValsMln_serial

Table containing the eigenvalues of the truncated response matrices from size `3*1` to `3*N`.
"""
tabEigValsMln_serial = Vector{Vector{Complex{Float64}}}([zeros(Float64,3*k) for k=1:N+1])

"""
    tabTruncMln!([args])

Fills the table `tabTruncMln` with the truncated matrices.
"""
function tabTruncMln!(tabTruncMln=tabTruncMln_serial, tabAln=tabAln_serial,
    tabBln=tabBln_serial, tabCln=tabCln_serial, tabDln=tabDln_serial,
    tabFln=tabFln_serial, tabGln=tabGln_serial, tabHln=tabHln_serial)
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

"""
    tabEigValsMln!([args])

Fills the table `tabTruncMln` with the truncated matrices.
"""
function tabEigValsMln!(tabTruncMln=tabTruncMln_serial,
    tabEigValsMln=tabEigValsMln_serial)
    for k=1:N+1
        tabEigValsMln[k] = eigvals(tabTruncMln[k])
    end
end

"""
    getPhysicalEigvals([args])

Computes the physical eigenvalues of the response matrix.
"""
function getPhysicalEigvals(tabEigValsMln=tabEigValsMln_serial)
    liEgv = zeros(Complex{Float64},3*(N+1))
    liEigvals = tabEigValsMln[N+1]
    liEgvCvCurrent = zeros(Complex{Float64},N)

    for iegv=1:3*(N+1) # take one of the eigenvalues of Mln(N+1)
        egv = liEigvals[iegv]
        for k=1:N # compare to those of its truncatures
            dist_min = -1
            egvSelect = 0
            for kegv=1:3*k # at a given truncature, take the closest one
                egvTrunc = tabEigValsMln[k][kegv]
                dist = abs(egv-egvTrunc)
                if ((dist_min<0) || (dist<dist_min))
                    egvSelect = egvTrunc
                    dist_min = dist
                end
            end
            if ((dist_min>=0) && (dist_min<error_dis)) # check if the closest one is within the accepted error margin
                liEgvCvCurrent[k] = egvSelect
            else # if not, put Inf to discard
                liEgvCvCurrent[k] = Inf
            end
        end

        nb_egv = 1 #(1 comes from initial separate egv at N+1)
        for k=1:N # count number of accepted truncature eigenvalues
            if (liEgvCvCurrent[k] != Inf)
                nb_egv += 1
            end
        end

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
                var_egv += abs2(liEgvCvCurrent[k]-moy_egv)
            end
        end
        var_egv = (var_egv+abs2(egv-moy_egv))/nb_egv

        if ((var_egv < error_var) && (nb_egv >= threshold_nb_egv)) # recover mean physical eigenvalue
            liEgv[iegv] = egv
        else
            liEgv[iegv] = Inf
        end

        for k=1:N
            liEgvCvCurrent[k] = 0.0
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

"""
    get_max(li)

Returns the tuple `(kEgv, egv)` where `egv` is the element of the list `li` with the highest
imaginary part and `kEgv` its position on the list.
"""
function get_max(li)
    len = length(li)
    max = -Inf
    egv = -Inf
    kEgv = 0
    for k=1:len
        if (imag(li[k])>max)
            egv = li[k]
            max = imag(egv)
            kEgv = k
        end
    end
    return kEgv, egv
end

######################################
# Clear table
######################################

"""
    tabTruncMln_clear!([args])

Clears the table containing the truncated response matrices.
"""
function tabTruncMln_clear!(tabTruncMln=tabTruncMln_serial)
    for k=1:N+1
        for l=1:k
            for n=1:k
                tabTruncMln[k][l,n] = 0.0
                tabTruncMln[k][l,k+n] = 0.0
                tabTruncMln[k][l,2*k+n] = 0.0
                tabTruncMln[k][k+l,n] = 0.0
                tabTruncMln[k][k+l,k+n] = 0.0
                tabTruncMln[k][k+l,2*k+n] = 0.0
                tabTruncMln[k][2*k+l,n] = 0.0
                tabTruncMln[k][2*k+l,k+n] = 0.0
                tabTruncMln[k][2*k+l,2*k+n] = 0.0
            end
        end
    end
end

"""
    tabEigValsMln_clear!([args])

Clears the table containing the eigenvalues.
"""
function tabEigValsMln_clear!(tabEigValsMln=tabEigValsMln_serial)
    for k=1:N+1
        tabEigValsMln[k] = zeros(Float64,3*k)
    end
end
