##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../sources/Main.jl") # Loading the main code

x = 0.0 # Mbulge/(Mdisk+Mbulb)
q = 1.0 # Mdisk/(Mdisk+Mdarkhalo)

println("Computing the eigenvalues...")

table_function_fill!(x,q)
matrix_fill!()
tabTruncMln!()
tabEigValsMln!(x,q)

@time physical_eig = getPhysicalEigvals()

########################################
namefile = "../data/Dump_Eigenvalues.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"tabEigVals",tabEigValsMln_serial[N+1])
    write(file,"tabEigValsPhys",physical_eig)
    close(file) # Closing the file
end

########################################
println("Saving data...")
writedump!(namefile) # Dumping the computed table

# clear the temp table
tabEigValsMln_clear!()
tabTruncMln_clear!()
matrix_clear!()
table_function_clear!()
