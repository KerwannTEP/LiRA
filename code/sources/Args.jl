##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--parallel"
    help = "Parallel computation: yes/no"
    arg_type = String
    default = "yes"
    "--eps"
    help = "Disk temperature"
    arg_type = Float64
    default = 0.1
    "--nbK"
    help = "Sampling number for integration"
    arg_type = Int64
    default = 5000
    "--N"
    help = "Number of basis elements"
    arg_type = Int64
    default = 170
    "--m"
    help = "Azimutal wavenumber"
    arg_type = Int64
    default = 2
end
parsed_args = parse_args(tabargs)

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
