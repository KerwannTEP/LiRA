##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../sources/Main.jl") # Loading the main code
########################################
# z=1.3
qminMeasure, qmaxMeasure = 0.01,1.0 # Range in q where the GR are computed
xminMeasure, xmaxMeasure = 0.0,0.4 # Range in x where the GR are computed

nbqMeasure = 10 # Number of q for which the GR are computed
nbxMeasure = 10 # Number of x for which the GR are computed
nbqxGrid = nbqMeasure*nbxMeasure
tabqMeasure = exp.(range(log(qminMeasure),length=nbqMeasure,log(qmaxMeasure)))
tabxMeasure = collect(range(xminMeasure,length=nbxMeasure,xmaxMeasure))


const tabqxGrid = zeros(Float64,2,nbqxGrid) # Location (q,x) of the grid points where the diffusion coefficients are computed
const tabGRGrid = zeros(Float64,nbqxGrid)
const tabRotGrid = zeros(Float64,nbqxGrid)


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


function tabGRRotGrid!()

    if (PARALLEL == "yes") # Computation is made in parallel
        println("parallel")
        Threads.@threads for iGrid=1:nbqxGrid
            q, x = tabqxGrid[1,iGrid], tabqxGrid[2,iGrid]

            tabOmega_parallel = zeros(Float64,nbK)
            tabKappaSqOverTwoOmega_parallel = zeros(Float64,nbK)

            tabAln_parallel = zeros(Float64,N+1,N+1)
            tabBln_parallel = zeros(Float64,N+1,N+1)
            tabCln_parallel = zeros(Float64,N+1,N+1)
            tabDln_parallel = zeros(Float64,N+1,N+1)
            tabFln_parallel = zeros(Float64,N+1,N+1)
            tabGln_parallel = zeros(Float64,N+1,N+1)
            tabHln_parallel = zeros(Float64,N+1,N+1)

            tabTruncMln_parallel = Vector{Array{Float64,2}}([zeros(Float64,3*k,3*k) for k=1:N+1])
            tabEigValsMln_parallel = Vector{Vector{Complex{Float64}}}([zeros(Float64,k) for k=1:N+1])

            re, im = growth_rate_with_rotation(x,q,tabOmega_parallel,tabKappaSqOverTwoOmega_parallel,
                            tabAln_parallel,tabBln_parallel,tabCln_parallel,tabDln_parallel,
                            tabFln_parallel,tabGln_parallel,tabHln_parallel,tabTruncMln_parallel,
                            tabEigValsMln_parallel)

            if (im > threshold_GR)
                tabRotGrid[iGrid], tabGRGrid[iGrid] = re, im
            else
                tabRotGrid[iGrid], tabGRGrid[iGrid] = 0.0, 0.0
            end
        end
    else # Computation is not made in parallel
        println("not parallel")

        for iGrid=1:nbqxGrid # Loop over the elements of the (a,j)-grid
            q, x = tabqxGrid[1,iGrid], tabqxGrid[2,iGrid] # Current (a,j) location

            re, im = growth_rate_with_rotation(x,q)

            if (im > threshold_GR)
                tabRotGrid[iGrid], tabGRGrid[iGrid] = re, im
            else
                tabRotGrid[iGrid], tabGRGrid[iGrid] = 0.0, 0.0
            end
        end
    end
end


########################################
namefile = "../data/Dump_Growth_Rate_Precession.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"a",a) # Characteristic length of disk
    write(file,"b",b) # Characteristic length of halo
    write(file,"c",c) # Characteristic length of bulge
    write(file,"eps",eps)
    write(file,"nbq",nbqMeasure) # Dumping the grid of (q,x)
    write(file,"nbx",nbxMeasure) # Dumping the grid of (q,x)
    write(file,"tabqx",tabqxGrid) # Dumping the grid of (q,x)
    write(file,"tabq",tabqMeasure) # Dumping the grid of q
    write(file,"tabx",tabxMeasure) # Dumping the grid of x
    write(file,"tabGR",tabGRGrid) # Dumping the total GR(q,x) for the current used value of lmax
    write(file,"tabRot",tabRotGrid) # Dumping the total Rot(q,x) for the current used value of lmax
    close(file) # Closing the file
end

########################################

tabqxGrid!()

println("Computing the maximum growth rate...")
@time tabGRRotGrid!()

########################################
println("Saving data...")
writedump!(namefile) # Dumping the computed table
