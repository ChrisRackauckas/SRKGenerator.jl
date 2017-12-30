using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
using ParameterizedFunctions, DiffEqMonteCarlo, Base.Test, OrdinaryDiffEq
using Plots; plotly()
using BenchmarkTools

################################################################################

prob = prob_sde_2Dlinear

reltols = 1.0./10.0.^(1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) + 2))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) + 2))
          Dict(:alg=>SRI())
          Dict(:alg=>SOSRI())
          Dict(:alg=>SOSRI2())
          ]
names = ["EM","RKMil","SRIW1","SOSRI","SOSRI2"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns=1000,names=names,error_estimate=:l2)

dts = 1./2.^(7:-1:4) #14->7 good plot
srand(100)
sim1 = test_convergence(dts,prob,SOSRI(),numMonte=10000)
sim2 = test_convergence(dts,prob,SOSRI2(),numMonte=10000)

p1 = plot(wp,title="Work-Precision Plot")
p2 = plot(dts,[sim1.errors[:l2],sim2.errors[:l2]],title="Convergence Tests",
          label=["SOSRI" "SOSRI2"], xscale = :log10, yscale = :log10, lw=3)
plot(p2,p1,layout=(2,1))
