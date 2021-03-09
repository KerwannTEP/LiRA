##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../../sources/ConstantStellarMass/Main.jl") # Loading the main code
########################################
xminMeasure, xmaxMeasure = 0.95,1.0 # Range in x where the GR are computed
qMeasure = 1.0 # Number of q for which the GR are computed
nbxMeasure = 100 # Number of x for which the GR are computed
tabxMeasure = collect(range(xminMeasure,length=nbxMeasure,xmaxMeasure))

const tabGR = zeros(Float64,nbxMeasure)

function tabGR!()
    for ix=1:nbxMeasure # Loop over the elements of the (a,j)-grid
        x = tabxMeasure[ix]
#        @time tabGR[ix,1], tabGR[ix,2]  = x, grow_rate(x,qMeasure)
        @time tabGR[ix]  = grow_rate(x,qMeasure)
    end
end



########################################

println("Computing the maximum growth rate...")
@time tabGR!()

########################################
println("Plotting the data...")
#p = Plots.plot(tabGR,xaxis="x",yaxis="GR")
p = Plots.plot(tabxMeasure,tabGR,xaxis="x",yaxis="GR")
Plots.savefig(p,"../../graphs/growth_rate_cut.png") # Saves the figure
Plots.display(p)
readline()
