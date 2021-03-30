"""
    growth_rate(x,q,[args])

Computes the growth rate of the system, i.e. the highest imaginary part of the physical eigenvales.
Empties the tables during the run.
"""
function growth_rate(x::Float64=0.0, q::Float64=1.0,
    tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabKappaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabKappaSqOverTwoOmega)
    tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    tabEigValsMln!(x,q,tabTruncMln,tabEigValsMln)
    physical_eig = getPhysicalEigvals(tabEigValsMln)
    kEgv, egv = get_max(physical_eig)

    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # return maximum growth rate
    return imag(egv)
end

"""
    growth_rate(x,q,[args])

Computes the physical eigenvalue with the highest imaginary part.
Empties the tables during the run.
"""
function growth_rate_with_rotation(x::Float64=0.0, q::Float64=1.0,
    tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabKappaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabKappaSqOverTwoOmega)
    tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    tabEigValsMln!(x,q,tabTruncMln,tabEigValsMln)
    physical_eig = getPhysicalEigvals(tabEigValsMln)
    kEgv, egv = get_max(physical_eig)

    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # return maximum growth rate
    return real(egv), imag(egv)
end

"""
    eigenmode_growth_rate(x,q,[args])

Computes the eigenvector corresponding to the physical eigenvalue with the highest imaginary part.
Empties the tables during the run.
"""
function eigenmode_growth_rate(x::Float64=0.0, q::Float64=1.0,
    tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabKappaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabKappaSqOverTwoOmega)
    tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    tabEigValsMln!(x,q,tabTruncMln,tabEigValsMln)
    physical_eig = getPhysicalEigvals(tabEigValsMln)
    kEgv, egv = get_max(physical_eig)

    kMin = 0
    minDist = Inf
    for i=1:3*(N+1)
        dist = abs(tabEigValsMln[N+1][i]-egv)
        if (dist < minDist)
            kMin = i
            minDist = dist
        end
    end

    # computes eigenmode of fastest growth rate
    egmode = tabEigVecsMln(x,q,tabTruncMln)[:,kMin]
    aln = egmode[1:N+1]
    bln = egmode[(N+1)+1:2*(N+1)]
    cln = egmode[2*(N+1)+1:3*(N+1)]

    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # return eigenmode of fastest growth rate
    return aln, bln, cln, egv
end
