include("../../sources/constant_stellar_mass/Main.jl") # Loading the main code

x = 1.0
q = 1.0


println("----------------------------------------------")

println("Computing the growth rate")
@time rate = grow_rate(x,q)
println("Growth rate  = ",rate)

println("----------------------------------------------")
