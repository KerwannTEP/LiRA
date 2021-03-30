"""
    physical_eigenvalues(x,q,[args])

Computes the physical eigenvalues of the system given the selection criteria.
Empties the tables during the run.
"""
function physical_eigenvalues(x::Float64=0.0, q::Float64=1.0,
    tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabKappaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabKappaSqOverTwoOmega)
    tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    tabEigValsMln!(tabTruncMln,tabEigValsMln)
    physical_eig = getPhysicalEigvals(tabEigValsMln)

    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # return maximum growth rate
    return physical_eig
end
