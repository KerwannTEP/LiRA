##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../../sources/constant_stellar_mass/Main.jl") # Loading the main code
########################################
qMeasure = 1.0 # q where the GR are computed
xminMeasure, xmaxMeasure = 0.5,1.0 # Range in x where the GR are computed
nbxMeasure = 20 # Number of x for which the GR are computed
tabxMeasure = collect(range(xminMeasure,length=nbxMeasure,xmaxMeasure))

const tabImagEigenMode = zeros(Float64,nbxMeasure)
const tabRealEigenMode = zeros(Float64,nbxMeasure)

function tabEigenMode!()

    if (PARALLEL == "yes") # Computation is made in parallel
        println("parallel")
        Threads.@threads for ix=1:nbxMeasure # Loop over the elements of the (a,j)-grid

            x = tabxMeasure[ix]

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

            tabRealEigenMode[ix], tabImagEigenMode[ix] = grow_rate_with_rotation(
                            x,qMeasure,tabOmega_parallel,tabKappaSqOverTwoOmega_parallel,
                            tabAln_parallel,tabBln_parallel,tabCln_parallel,tabDln_parallel,
                            tabFln_parallel,tabGln_parallel,tabHln_parallel,tabTruncMln_parallel,
                            tabEigValsMln_parallel)

        end
    else # Computation is not made in parallel
        println("not parallel")

        for ix=1:nbxMeasure # Loop over the elements of the (a,j)-grid
            x = tabxMeasure[ix]
            tabRealEigenMode[ix], tabImagEigenMode[ix] = grow_rate_with_rotation(x,qMeasure)
        end
    end
end


########################################
namefile = "../../data/Dump_Eigenmode_eps0_1_q_1.0.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"q",qMeasure) # Dumping the grid of (q,x)
    write(file,"nbx",nbxMeasure) # Dumping the grid of (q,x)
    write(file,"tabx",tabxMeasure) # Dumping the grid of x
    write(file,"tabRealEigenMode",tabRealEigenMode) # Dumping the total GR(q,x) for the current used value of lmax
    write(file,"tabImagEigenMode",tabImagEigenMode) # Dumping the total Rot(q,x) for the current used value of lmax
    close(file) # Closing the file
end

########################################


println("Computing the eigenmodes...")
@time tabEigenMode!()

########################################
println("Saving data...")
writedump!(namefile) # Dumping the computed table
