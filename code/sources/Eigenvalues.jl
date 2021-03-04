using Plots

function physical_eigenvalues(x::Float64=1.0, q::Float64=1.0)
    table_function_fill!(x,q)
    matrix_fill!()
    tabTruncMln!()
    tabEigValsMln!(x,q)
    physical_eig = getPhysicalEigvals()

    # clear the temp table
    tabEigValsMln_clear!()
    tabTruncMln_clear!()
    matrix_clear!()
    table_function_clear!()

    # return maximum growth rate
    return physical_eig
end

function plot_eigenvalues(x::Float64=1.0,q::Float64=1.0)
    table_function_fill!(x,q)
    matrix_fill!()
    tabTruncMln!()
    tabEigValsMln!(x,q)
    eig = tabEigValsMln[N+1]

    loci_egv_real = zeros(Float64,3*(N+1)) #(real part, imag part)
    loci_egv_imag = zeros(Float64,3*(N+1))

    for k=1:3*(N+1)
        loci_egv_real[k], loci_egv_imag[k]   = real(eig[k]) ,imag(eig[k])
    end
    # clear the temp table
    tabEigValsMln_clear!()
    tabTruncMln_clear!()
    matrix_clear!()
    table_function_clear!()

    # plot
    p = Plots.scatter(loci_egv_real, [loci_egv_imag])
    Plots.display(p)
    readline()
#    return loci_egv
end

function plot_phys_eigenvalues(x::Float64=1.0,q::Float64=1.0)
    table_function_fill!(x,q)
    matrix_fill!()
    tabTruncMln!()
    tabEigValsMln!(x,q)
    physical_eig = getPhysicalEigvals()

    nb_egv = length(physical_eig)
    loci_egv_real = zeros(Float64,3*(N+1)) #(real part, imag part)
    loci_egv_imag = zeros(Float64,3*(N+1))

    for k=1:nb_egv
        loci_egv_real[k], loci_egv_imag[k] = real(physical_eig[k]), imag(physical_eig[k])
    end
    # clear the temp table
    tabEigValsMln_clear!()
    tabTruncMln_clear!()
    matrix_clear!()
    table_function_clear!()

    # plot
    p = Plots.scatter(loci_egv_real, [loci_egv_imag])
    Plots.display(p)
    readline()
#    return loci_egv
end
