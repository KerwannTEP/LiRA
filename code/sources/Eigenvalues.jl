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
