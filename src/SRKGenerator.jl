__precompile__()

module SRKGenerator

using NLopt, CUDArt

const SOFT_C = true

include("getCoefs.jl")
include("main.jl")
include("optim_functions.jl")
include("printing.jl")
include("translate.jl")

export srk_optimize

end # module
