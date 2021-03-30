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

"""
    plot_eigenvalues(x,q,[args])

Plots the eigenvalues of the truncated response matrix on a Nyquist diagramm.
Empties the tables during the run.
"""
function plot_eigenvalues(x::Float64=0.0,q::Float64=1.0,
    tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabKappaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabKappaSqOverTwoOmega)
    tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    tabEigValsMln!(tabTruncMln,tabEigValsMln)
    eig = tabEigValsMln[N+1]

    loci_egv_real = zeros(Float64,3*(N+1))
    loci_egv_imag = zeros(Float64,3*(N+1))

    for k=1:3*(N+1)
        loci_egv_real[k], loci_egv_imag[k]   = real(eig[k]) ,imag(eig[k])
    end
    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # plot
    p = Plots.scatter(loci_egv_real, [loci_egv_imag]
    ,xlim=(-5,10),ylim=(-2.0,2.0))

    Plots.display(p)
    readline()
end

"""
    plot_phys_eigenvalues(x,q,[args])

Plots the physical eigenvalues of the response matrix on a Nyquist diagramm.
Empties the tables during the run.
"""
function plot_phys_eigenvalues(x::Float64=0.0,q::Float64=1.0,
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

    nb_egv = length(physical_eig)
    loci_egv_real = zeros(Float64,3*(N+1))
    loci_egv_imag = zeros(Float64,3*(N+1))

    for k=1:nb_egv
        loci_egv_real[k], loci_egv_imag[k] = real(physical_eig[k]), imag(physical_eig[k])
    end
    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # plot
    p = Plots.scatter(loci_egv_real,[loci_egv_imag])
    Plots.display(p)
    readline()
end
