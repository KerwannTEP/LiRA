function grow_rate(x::Float64=1.0, q::Float64=1.0,
    tabOmega=tabOmega_serial, tabAlphaSqOverTwoOmega=tabAlphaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabAlphaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabAlphaSqOverTwoOmega)
    tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    tabEigValsMln!(x,q,tabTruncMln,tabEigValsMln)
    physical_eig = getPhysicalEigvals(tabEigValsMln)
    kEgv, egv = get_max(physical_eig)

    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabAlphaSqOverTwoOmega)

    # return maximum growth rate
    return imag(egv)
end

function eigenmode_grow_rate(x::Float64=1.0, q::Float64=1.0,
    tabOmega=tabOmega_serial, tabAlphaSqOverTwoOmega=tabAlphaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabAlphaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabAlphaSqOverTwoOmega)
    tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    tabEigValsMln!(x,q,tabTruncMln,tabEigValsMln)
    physical_eig = getPhysicalEigvals(tabEigValsMln)
    kEgv, egv = get_max(physical_eig)

    # computes eigenmode of fastest growth rate
    egmode = tabEigVecsMln(x,q,tabTruncMln)[:,kEgv]
    aln = egmode[1:N+1]
    bln = egmode[(N+1)+1:2*(N+1)]
    cln = egmode[2*(N+1)+1:3*(N+1)]

    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabAlphaSqOverTwoOmega)

    # return eigenmode of fastest growth rate
    return aln, bln, cln
end
