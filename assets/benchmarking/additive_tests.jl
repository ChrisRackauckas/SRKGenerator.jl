using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
using ParameterizedFunctions, DiffEqMonteCarlo, Base.Test, OrdinaryDiffEq
using Plots; gr()
using BenchmarkTools

################################################################################

## Additive

prob = prob_sde_additive

reltols = 1.0./10.0.^(-2:4)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) - 2))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) - 2))
          Dict(:alg=>SRA1())
          Dict(:alg=>SRA2())
          Dict(:alg=>SRA3())
          Dict(:alg=>SOSRA())
          Dict(:alg=>SOSRA2())
          Dict(:alg=>RackKenCarp())
          ]
names = ["EM","RKMil","SRA1","SRA2","SRA3","SOSRA","SOSRA2","SKenCarp"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns=10000,names=names,error_estimate=:l2)

dts = 1./2.^(10:-1:2) #14->7 good plot
sim4 = test_convergence(dts,prob,SOSRA(),numMonte=Int(1e1))
sim5 = test_convergence(dts,prob,SOSRA2(),numMonte=Int(1e1))
sim6 = test_convergence(dts,prob,RackKenCarp(),numMonte=Int(1e1))

p1 = plot(wp,title="Work-Precision Plot")
p2 = plot(dts,[sim4.errors[:l2],sim5.errors[:l2],sim6.errors[:l2]],
          title="Convergence Tests",label=["SOSRA" "SOSRA2" "SKenCarp"],
          yscale=:log10,xscale=:log10,lw=3)
plot(p2,p1,layout=(2,1))
