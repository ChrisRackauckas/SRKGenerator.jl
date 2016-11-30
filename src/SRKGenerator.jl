module SRKGenerator

using NLopt, CUDArt

include("getCoefs.jl")
include("main.jl")
include("optim_functions.jl")
include("printing.jl")

export SRKoptimize

end # module
