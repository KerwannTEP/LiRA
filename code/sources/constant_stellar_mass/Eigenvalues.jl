using Plots

function physical_eigenvalues(x::Float64=0.0, q::Float64=1.0,
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

    # clear the temp table
    tabEigValsMln_clear!(tabEigValsMln)
    tabTruncMln_clear!(tabTruncMln)
    matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # return maximum growth rate
    return physical_eig
end

# Returns a table of the nB = 11 eigenvalues of the B-branch (see Aoki 1979) for x=0 and q=1
# function initial_physical_eigenvalues(
#         tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial,
#         tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
#         tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
#         tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
#         tabEigValsMln=tabEigValsMln_serial)
#
#     x = 0.0
#     q = 1.0
#     nB = 11
#     table_function_fill!(x,q,tabOmega,tabKappaSqOverTwoOmega)
#     matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabKappaSqOverTwoOmega)
#     tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
#     tabEigValsMln!(x,q,tabTruncMln,tabEigValsMln)
#     physical_eig = getPhysicalEigvals(tabEigValsMln)
#
#     tabBmodes = zeros(Complex{Float64},nB)
#     n_modes = length(physical_eig)
#
#     mode = 0.0
#
#     # fill max eigenvalue
#
#     kEgv, egv = get_max(physical_eig)
#     tabBmodes[1] = egv
#
#     for i=2:nB
#         mode = 0.0
#         # Fill max eigenvalues
#         for jmode=1:n_modes
#             if ((imag(physical_eig[jmode])>imag(mode)) && (imag(physical_eig[jmode])<imag(tabBmodes[i-1])))
#                 mode = physical_eig[jmode]
#             end
#
#         # we have the next highest growth rate
#         tabBmodes[i] = mode
#         end
#     end
#
#
#
#
#
#     # clear the temp table
#     tabEigValsMln_clear!(tabEigValsMln)
#     tabTruncMln_clear!(tabTruncMln)
#     matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
#     table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)
#
#     # return maximum growth rate
#     return tabBmodes
# end

function plot_eigenvalues(x::Float64=0.0,q::Float64=1.0,
    tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial,
        tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
        tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
        tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
        tabEigValsMln=tabEigValsMln_serial)
    table_function_fill!(x,q,tabOmega,tabKappaSqOverTwoOmega)
    matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabKappaSqOverTwoOmega)
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
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # plot
    p = Plots.scatter(loci_egv_real, [loci_egv_imag]
    ,xlim=(-5,10),ylim=(-2.0,2.0))

    Plots.display(p)
    readline()
#    return loci_egv
end

function plot_phys_eigenvalues(x::Float64=0.0,q::Float64=1.0,
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
    table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)

    # plot
    p = Plots.scatter(loci_egv_real, [loci_egv_imag]
    )#,xlim=(-5,10),ylim=(-0.25,0.25))
    Plots.display(p)
    readline()
#    return loci_egv
end

# function plot_phys_eigenvalues_adv(x::Float64=0.0,q::Float64=1.0,
#     tabOmega=tabOmega_serial, tabKappaSqOverTwoOmega=tabKappaSqOverTwoOmega_serial,
#         tabAln=tabAln_serial, tabBln=tabBln_serial, tabCln=tabCln_serial,
#         tabDln=tabDln_serial, tabFln=tabFln_serial, tabGln=tabGln_serial,
#         tabHln=tabHln_serial, tabTruncMln=tabTruncMln_serial,
#         tabEigValsMln=tabEigValsMln_serial, tabEigVecsMln=tabEigVecsMln_serial)
#     table_function_fill!(x,q,tabOmega,tabKappaSqOverTwoOmega)
#     matrix_fill!(x,q,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln,tabOmega,tabKappaSqOverTwoOmega)
#     tabTruncMln!(tabTruncMln,tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
#     tabEigValsMln!(x,q,tabTruncMln,tabEigValsMln)
#     tabEigVecsMln!(x,q,tabTruncMln,tabEigVecsMln)
#     physical_eig = getPhysicalEigvalsAdvanced(tabEigValsMln,tabEigVecsMln)
#
#     nb_egv = length(physical_eig)
#     loci_egv_real = zeros(Float64,3*(N+1)) #(real part, imag part)
#     loci_egv_imag = zeros(Float64,3*(N+1))
#
#     for k=1:nb_egv
#         loci_egv_real[k], loci_egv_imag[k] = real(physical_eig[k]), imag(physical_eig[k])
#     end
#     # clear the temp table
#     tabEigVecsMln_clear!(tabEigVecsMln)
#     tabEigValsMln_clear!(tabEigValsMln)
#     tabTruncMln_clear!(tabTruncMln)
#     matrix_clear!(tabAln,tabBln,tabCln,tabDln,tabFln,tabGln,tabHln)
#     table_function_clear!(tabOmega,tabKappaSqOverTwoOmega)
#
#     # plot
#     p = Plots.scatter(loci_egv_real, [loci_egv_imag]
#     )#,xlim=(-5,10),ylim=(-0.25,0.25))
#     Plots.display(p)
#     readline()
# #    return loci_egv
# end
