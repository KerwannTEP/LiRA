include("Packages.jl") # Import the packages
include("Args.jl") # Parsing the command-line

##################################################
# General parameters
##################################################
"""
    PARALLEL

Determining if the code is run in parallel.
Is 'yes' by default.
"""
const PARALLEL = parsed_args["parallel"]
if ((PARALLEL != "yes") && (PARALLEL != "no"))
    error("ERROR: UNKNOWN PARALLEL") # Unknown parallel procedure
end


"""
    eps

Temperature of the disk, defined as the ratio of the internal
energy over the absolut total energy.

A low temperature means that pressure predominates, while a high
temperature means that self-gravity predominates.

Is '0.1' by default.
"""
const eps = parsed_args["eps"]

"""
    N

Number of basis elements in the linear expansion.

Is '100' by default.
"""
const N = parsed_args["N"]

"""
    nbK

Sampling number for numerical integration.

Is '5000' by default.
"""
const nbK = parsed_args["nbK"]

"""
    m

Azimutal wavenumber in the Fourier expansion in angles.

Is '2' by default (mode with fastest growth rates).
"""
const m = parsed_args["m"]

##################################################

include("Constants.jl") # Physical constants
include("Mean.jl") # Useful functions of the system
include("TableValues.jl") # Compute the matrix elements
include("ResponseMatrix.jl") # Compute the linear response matrix
include("Eigenvalues.jl") # Computes the system eigenvalues
include("GrowthRate.jl") # Computes the growth rate of the system
