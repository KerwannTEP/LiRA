##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../sources/Main.jl") # Loading the main code

x = 1.0 # Mdisk/(Mdisk+Mbulb)
q = 1.0 # Mdisk/(Mdisk+Mdarkhalo)

println("Computing the eigenvalues...")

aln, bln, cln = eigenmode_grow_rate(x,q)

########################################
namefile = "../data/Dump_Eigenmode.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"aln",aln)
    write(file,"bln",bln)
    write(file,"cln",cln)
    close(file) # Closing the file
end

########################################
println("Saving data...")
writedump!(namefile) # Dumping the computed table
