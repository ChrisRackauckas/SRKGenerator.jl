using StochasticDiffEq, DiffEqProblemLibrary
prob = oval2ModelExample(largeFluctuations=true,useBigs=false)

# Most are found at DiffEqBenchmarks.jl

srand(250)
@time sol = solve(prob,ImplicitEM(),dt=1/10000,progress=true)
