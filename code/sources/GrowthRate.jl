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
    tabEigValsMln!(tabTruncMln,tabEigValsMln)
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
    growth_rate_with_rotation(x,q,[args])

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
    tabEigValsMln!(tabTruncMln,tabEigValsMln)
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
