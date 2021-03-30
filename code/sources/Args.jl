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
