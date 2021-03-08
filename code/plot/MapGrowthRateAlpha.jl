##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../sources/Main.jl") # Loading the main code
########################################
aminMeasure, amaxMeasure = 0.0,10.0 # Range in q where the GR are computed
xminMeasure, xmaxMeasure = 0.1,1.0 # Range in x where the GR are computed
nbaMeasure = 200 # Number of q for which the GR are computed
nbxMeasure = 200 # Number of x for which the GR are computed
nbaxGrid = nbaMeasure*nbxMeasure # Number of (a,j) for which the Djj are computed
tabaMeasure = collect(range(aminMeasure,length=nbaMeasure,amaxMeasure))
tabxMeasure = collect(range(xminMeasure,length=nbxMeasure,xmaxMeasure))

const tabaxGrid = zeros(Float64,2,nbaxGrid) # Location (q,x) of the grid points where the diffusion coefficients are computed
const tabGRGrid = zeros(Float64,nbaxGrid)

PARALLEL = "yes"


function tabaxGrid!()
    index = 1
    for ia=1:nbaMeasure
    aMeasure = tabaMeasure[ia]
        for ix=1:nbxMeasure
            xMeasure = tabxMeasure[ix]
            tabaxGrid[1,index], tabaxGrid[2,index] = aMeasure, xMeasure
            index += 1
        end
    end
end

function init_tabGRGrid!()
    for iGrid=1:nbaxGrid
        tabGRGrid[iGrid] = 0.0
    end
end

function tabGRGrid!()
    init_tabGRGrid!() # Making sure that the grid is initially set to 0

    if (PARALLEL == "yes") # Computation is made in parallel
        println("parallel")
        Threads.@threads for iGrid=1:nbaxGrid # Loop over the elements of the (a,j)-grid
            a, x = tabaxGrid[1,iGrid], tabaxGrid[2,iGrid] # Current (a,j) location
            q = _q(a)
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

        for iGrid=1:nbaxGrid # Loop over the elements of the (a,j)-grid
            a, x = tabaxGrid[1,iGrid], tabaxGrid[2,iGrid] # Current (a,j) location
            q = _q(a)

            tabGRGrid[iGrid] = grow_rate(x,q)
        end
    end
end


########################################
namefile = "../data/Dump_Growth_Rate_alpha_eps0_2.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"nba",nbaMeasure) # Dumping the grid of (q,x)
    write(file,"nbx",nbxMeasure) # Dumping the grid of (q,x)
    write(file,"tabax",tabaxGrid) # Dumping the grid of (q,x)
    write(file,"taba",tabaMeasure) # Dumping the grid of q
    write(file,"tabx",tabxMeasure) # Dumping the grid of x
    write(file,"tabGR",tabGRGrid) # Dumping the total GR(q,x) for the current used value of lmax
    close(file) # Closing the file
end

########################################

tabaxGrid!()

println("Computing the maximum growth rate...")
@time tabGRGrid!()

########################################
println("Saving data...")
writedump!(namefile) # Dumping the computed table
