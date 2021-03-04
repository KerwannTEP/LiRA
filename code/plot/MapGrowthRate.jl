##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../sources/Main.jl") # Loading the main code
########################################
qminMeasure, qmaxMeasure = 0.1,1.0 # Range in q where the GR are computed
xminMeasure, xmaxMeasure = 0.8,1.0 # Range in x where the GR are computed
nbqMeasure = 100 # Number of q for which the GR are computed
nbxMeasure = 100 # Number of x for which the GR are computed
nbqxGrid = nbqMeasure*nbxMeasure # Number of (a,j) for which the Djj are computed
tabqMeasure = collect(range(qminMeasure,length=nbqMeasure,qmaxMeasure))
tabxMeasure = collect(range(xminMeasure,length=nbxMeasure,xmaxMeasure))

const tabqxGrid = zeros(Float64,2,nbqxGrid) # Location (q,x) of the grid points where the diffusion coefficients are computed
const tabGRGrid = zeros(Float64,nbqxGrid)


function tabqxGrid!()
    index = 1
    for iq=1:nbqMeasure
    qMeasure = tabqMeasure[iq]
        for ix=1:nbxMeasure
            xMeasure = tabxMeasure[ix]
            tabqxGrid[1,index], tabqxGrid[2,index] = qMeasure, xMeasure
            index += 1
        end
    end
end

function init_tabGRGrid!()
    for iGrid=1:nbqxGrid
        tabGRGrid[iGrid] = 0.0
    end
end

function tabGRGrid!()
    init_tabGRGrid!() # Making sure that the grid is initially set to 0
    for iGrid=1:nbqxGrid # Loop over the elements of the (a,j)-grid
        q, x = tabqxGrid[1,iGrid], tabqxGrid[2,iGrid] # Current (a,j) location
        tabGRGrid[iGrid] = grow_rate(x,q)
    end
end


########################################
namefile = "../data/Dump_Growth_Rate_zoom.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"nbq",nbqMeasure) # Dumping the grid of (q,x)
    write(file,"nbx",nbxMeasure) # Dumping the grid of (q,x)
    write(file,"tabqx",tabqxGrid) # Dumping the grid of (q,x)
    write(file,"tabq",tabqMeasure) # Dumping the grid of q
    write(file,"tabx",tabxMeasure) # Dumping the grid of x
    write(file,"tabGR",tabGRGrid) # Dumping the total GR(q,x) for the current used value of lmax
    close(file) # Closing the file
end

########################################

tabqxGrid!()

println("Computing the maximum growth rate...")
@time tabGRGrid!()

########################################
println("Saving data...")
writedump!(namefile) # Dumping the computed table
