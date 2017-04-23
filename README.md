# SRKGenerator

This is the repository for the SRKGenerator of high strong order
Stochastic Runge-Kutta methods.

## Installation

To install this package, use

```julia
Pkg.clone("https://github.com/ChrisRackauckas/SRKGenerator.jl")
```

You may also need to run

```julia
Pkg.add("NLopt")
Pkg.add("CUDArt")
```

The CUDA runtime needs to be installed separately, and the CUDA kernal should
be compiled to a .ptx and added to the `/deps` folder.

## Running The Generator

The generator is setup on the `srk_optimize` the function. However, it's easiest
to run it as a test. See `test/runtests.jl` to modify the script, and use

```julia
Pkg.test("SRKGenerator")
```

to run the test script. See `src/main.jl` for the generator function and its
associated options. This method will make it re-precompile, so it's usually
better to run the script directly:

```julia
include(joinpath(Pkg.dir("SRKGenerator"),"test","runtests.jl"))
```

A command line utility is setup so that one can do:

```julia
julia cudaCores imin jmax
```

For example, for a highly deterministically stable method:

```julia
julia runtests.jl 2496 -12 1
```

For batch scripts, see the `/scripts` directory.

## Generator Options

```julia
function srk_optimize(alg,dx,pop_size,imin,imax,jmin,jmax;
                    NLoptRandSeed = 0,parameter_minmax=5,max_eval=Int(1e8),
                    initCon = ones(44),tol = 1e-2,ftol = 1e-15,tol2 = 1e-5,
                    counterSteps=Int(1e5),counterSteps2=Int(1e6),
                    initStepSize=[],gpuEnabled=true,ptx_str  = "integration.ptx",
                    cudaCores = 1664,initStepSize2=1e-6,outfile="",constrain_c = true)
```

The main options for the generator are:

- `alg`:  The algorithm which NLopt should use. This algorithm needs to be derivative-free
  and allow equality constraints. Possible choices are thus `:LN_AUGLAG_EQ`,
  `:LN_COBYLA`, and `:GN_ISRES`.
- `dx`: The `dx` in the integration. Smaller values are more precise but take
  longer to perform calculations
- `pop_size`: Controls the population size for `:GN_ISRES`
- `imin`,`imax`,`jmin`,`jmax`: Controls the integration boundaries.
- `NLoptRandSeed`: The random number seed for NLopt.
- `cudaCores`: The number of CUDA cores on your device.
- `ptx_str`: The path to the .ptx CUDA kernal.
