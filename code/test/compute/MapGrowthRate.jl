##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../../sources/constant_stellar_mass/Main.jl") # Loading the main code
########################################
qminMeasure, qmaxMeasure = 1.0,10.0 # Range in q where the GR are computed
xminMeasure, xmaxMeasure = 0.001,1.0 # Range in x where the GR are computed
nbqMeasure = 20 # Number of q for which the GR are computed
nbxMeasure = 20 # Number of x for which the GR are computed
nbqxGrid = nbqMeasure*nbxMeasure # Number of (a,j) for which the Djj are computed
tabqMeasure = collect(range(qminMeasure,length=nbqMeasure,qmaxMeasure))
tabxMeasure = collect(range(xminMeasure,length=nbxMeasure,xmaxMeasure))

const tabqxGrid = zeros(Float64,2,nbqxGrid) # Location (q,x) of the grid points where the diffusion coefficients are computed
const tabGRGrid = zeros(Float64,nbqxGrid)

PARALLEL = "yes"


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

    if (PARALLEL == "yes") # Computation is made in parallel
        println("parallel")
        Threads.@threads for iGrid=1:nbqxGrid # Loop over the elements of the (a,j)-grid
            q, x = tabqxGrid[1,iGrid], tabqxGrid[2,iGrid] # Current (a,j) location

            tabOmega_parallel = zeros(Float64,nbK)
            tabAlphaSqOverTwoOmega_parallel = zeros(Float64,nbK)

            tabAln_parallel = zeros(Float64,N+1,N+1)
            tabBln_parallel = zeros(Float64,N+1,N+1)
            tabCln_parallel = zeros(Float64,N+1,N+1)
            tabDln_parallel = zeros(Float64,N+1,N+1)
            tabFln_parallel = zeros(Float64,N+1,N+1)
            tabGln_parallel = zeros(Float64,N+1,N+1)
            tabHln_parallel = zeros(Float64,N+1,N+1)

            tabTruncMln_parallel = Vector{Array{Float64,2}}([zeros(Float64,3*k,3*k) for k=1:N+1])
            tabEigValsMln_parallel = Vector{Vector{Complex{Float64}}}([zeros(Float64,k) for k=1:N+1])

            tabGRGrid[iGrid] = grow_rate(x,q,tabOmega_parallel,tabAlphaSqOverTwoOmega_parallel,
                            tabAln_parallel,tabBln_parallel,tabCln_parallel,tabDln_parallel,
                            tabFln_parallel,tabGln_parallel,tabHln_parallel,tabTruncMln_parallel,
                            tabEigValsMln_parallel)

        end
    else # Computation is not made in parallel
        println("not parallel")

        for iGrid=1:nbqxGrid # Loop over the elements of the (a,j)-grid
            q, x = tabqxGrid[1,iGrid], tabqxGrid[2,iGrid] # Current (a,j) location

            tabGRGrid[iGrid] = grow_rate(x,q)
        end
    end
end


########################################
namefile = "../../data/Dump_Growth_Rate_eps0_1.hf5"
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
