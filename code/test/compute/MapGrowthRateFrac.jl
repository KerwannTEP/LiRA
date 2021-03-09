##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../../sources/constant_stellar_mass/Main.jl") # Loading the main code
########################################
fDHminMeasure, fDHmaxMeasure = 0.0,10.0 # Range in q where the GR are computed
xminMeasure, xmaxMeasure = 0.001,1.0 # Range in x where the GR are computed
nbfDHMeasure = 100 # Number of q for which the GR are computed
nbxMeasure = 100 # Number of x for which the GR are computed
nbfracGrid = nbfDHMeasure*nbxMeasure # Number of (a,j) for which the Djj are computed
tabfDHMeasure = collect(range(fDHminMeasure,length=nbfDHMeasure,fDHmaxMeasure))
tabxMeasure = collect(range(xminMeasure,length=nbxMeasure,xmaxMeasure))

const tabfracGrid = zeros(Float64,2,nbfracGrid) # Location (q,x) of the grid points where the diffusion coefficients are computed
const tabGRGrid = zeros(Float64,nbfracGrid)

PARALLEL = "yes"


function tabfracGrid!()
    index = 1
    for ifDH=1:nbfDHMeasure
    fDHMeasure = tabfDHMeasure[ifDH]
        for ix=1:nbxMeasure
            xMeasure = tabxMeasure[ix]
            tabfracGrid[1,index], tabfracGrid[2,index] = fDHMeasure, xMeasure
            index += 1
        end
    end
end

function init_tabGRGrid!()
    for iGrid=1:nbfracGrid
        tabGRGrid[iGrid] = 0.0
    end
end

function tabGRGrid!()
    init_tabGRGrid!() # Making sure that the grid is initially set to 0

    if (PARALLEL == "yes") # Computation is made in parallel
        println("parallel")
        Threads.@threads for iGrid=1:nbfracGrid # Loop over the elements of the (a,j)-grid
            fDH, x = tabfracGrid[1,iGrid], tabfracGrid[2,iGrid] # Current (a,j) location
            q = _q(x,fDH)

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

        for iGrid=1:nbfracGrid # Loop over the elements of the (a,j)-grid
            fDH, x = tabfracGrid[1,iGrid], tabfracGrid[2,iGrid] # Current (a,j) location
            q = _q(x,fDH)

            tabGRGrid[iGrid] = grow_rate(x,q)
        end
    end
end


########################################
namefile = "../../data/Dump_Growth_Rate_frac_eps0_1.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"nbfDH",nbfDHMeasure) # Dumping the grid of (q,x)
    write(file,"nbx",nbxMeasure) # Dumping the grid of (q,x)
    write(file,"tabfrac",tabfracGrid) # Dumping the grid of (q,x)
    write(file,"tabfDH",tabfDHMeasure) # Dumping the grid of q
    write(file,"tabx",tabxMeasure) # Dumping the grid of x
    write(file,"tabGR",tabGRGrid) # Dumping the total GR(q,x) for the current used value of lmax
    close(file) # Closing the file
end

########################################

tabfracGrid!()

println("Computing the maximum growth rate...")
@time tabGRGrid!()

########################################
println("Saving data...")
writedump!(namefile) # Dumping the computed table
