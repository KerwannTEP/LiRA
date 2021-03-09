using Plots

function physical_eigenvalues(x::Float64=1.0, q::Float64=1.0,
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

    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabAlphaSqOverTwoOmega)

    # return maximum growth rate
    return physical_eig
end

function plot_eigenvalues(x::Float64=1.0,q::Float64=1.0,
    tabOmega=tabOmega_serial, tabAlphaSqOverTwoOmega=tabAlphaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabAlphaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabAlphaSqOverTwoOmega)
    tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    tabEigValsMln!(x,q,tabTruncMln,tabEigValsMln)
    eig = tabEigValsMln[N+1]

    loci_egv_real = zeros(Float64,3*(N+1)) #(real part, imag part)
    loci_egv_imag = zeros(Float64,3*(N+1))

    for k=1:3*(N+1)
        loci_egv_real[k], loci_egv_imag[k]   = real(eig[k]) ,imag(eig[k])
    end
    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabAlphaSqOverTwoOmega)

    # plot
    p = Plots.scatter(loci_egv_real, [loci_egv_imag])
    Plots.display(p)
    readline()
#    return loci_egv
end

function plot_phys_eigenvalues(x::Float64=1.0,q::Float64=1.0,
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

    nb_egv = length(physical_eig)
    loci_egv_real = zeros(Float64,3*(N+1)) #(real part, imag part)
    loci_egv_imag = zeros(Float64,3*(N+1))

    for k=1:nb_egv
        loci_egv_real[k], loci_egv_imag[k] = real(physical_eig[k]), imag(physical_eig[k])
    end
    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabAlphaSqOverTwoOmega)

    # plot
    p = Plots.scatter(loci_egv_real, [loci_egv_imag])
    Plots.display(p)
    readline()
#    return loci_egv
end
