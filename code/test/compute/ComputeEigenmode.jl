##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../../sources/constant_stellar_mass/Main.jl") # Loading the main code

x = 0.0 # Mbulge/(Mdisk+Mbulb)
q = 1.0 # Mdisk/(Mdisk+Mdarkhalo)

println("Computing the eigenvalues...")

aln, bln, cln, egv = eigenmode_grow_rate(x,q)
println(egv)

########################################
namefile = "../../data/Dump_Eigenmode.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"a",a) # Characteristic length of disk
    write(file,"c",c) # Characteristic length of bulge
    write(file,"x",x)
    write(file,"q",q) 
    write(file,"aln",aln)
    write(file,"bln",bln)
    write(file,"cln",cln)
    write(file,"m",m)
    write(file,"egv",egv)
    write(file,"eps",eps)
    close(file) # Closing the file
end

########################################
println("Saving data...")
writedump!(namefile) # Dumping the computed table
