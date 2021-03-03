using HDF5 # To have access to .hf5 files
using Plots # To be able to plot data

########################################
# Function that read the .hf5 file
########################################

namefile = "../data/Dump_Growth_Rate_high_q.hf5"

"""
    openData(namefile)

Recovers the data from the .hf5 file `namefile` and returns it as a 4-uple.
In the order nbj, aMeasure, tabj, tabDRRjj.
"""
function openData(namefile)
    file = h5open(namefile, "r")
    nbq = read(file,"nbq")
    nbx = read(file,"nbx")
    tabqx = read(file,"tabqx")
    tabq = read(file,"tabq")
    tabx = read(file,"tabx")
    tabGR = read(file,"tabGR")
    close(file)
    return nbq, nbx, tabqx, tabq, tabx, tabGR
end

########################################
# Getting the data
########################################
println("Recovering plot data...")
nbq, nbx, tabqx, tabq, tabx, tabGR = openData(namefile)

#println(tabq)
#println(tabx)
#println(tabGR)

########################################
# Plotting the data
########################################
println("Plotting the data...")
p = Plots.contourf(tabq,tabx,tabGR,xaxis="x",yaxis="q")
Plots.savefig(p,"../graphs/growth_rate_high_q.png") # Saves the figure
Plots.display(p)
readline()
