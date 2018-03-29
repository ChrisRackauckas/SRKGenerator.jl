using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
using ParameterizedFunctions, DiffEqMonteCarlo, Base.Test, OrdinaryDiffEq
using Plots; pyplot()
using BenchmarkTools

################################################################################

## Additive

prob = prob_sde_additive

dts = 1./2.^(10:-1:2) #14->7 good plot
sim1 = test_convergence(dts,prob,SOSRA(),numMonte=Int(1e3))
sim2 = test_convergence(dts,prob,SOSRA2(),numMonte=Int(1e3))
sim3 = test_convergence(dts,prob,RackKenCarp(),numMonte=Int(1e3))

p1 = plot(dts,[sim1.errors[:l2],sim2.errors[:l2],sim3.errors[:l2]],
          title="Additive Convergence Tests",label=["SOSRA" "SOSRA2" "SKenCarp"],
          xtickfont = font(16, "LM Roman"),titlefont = font(20, "LM Roman"),
          ytickfont = font(16, "LM Roman"),guidefont = font(18, "LM Roman"),
          legendfont = font(14, "LM Roman"),
          ylabel = "Error", xlabel="dt",
          yscale=:log10,xscale=:log10,lw=3)

dts = 1./2.^(10:-1:2) #14->7 good plot
ref = dts.^2 * 1e-4
plot!(p1,dts,ref,linestyle=:dash,lw=3,label="2nd Order Reference")

prob = prob_sde_2Dlinear

dts = 1./2.^(7:-1:4) #14->7 good plot
srand(100)
sim1 = test_convergence(dts,prob,SOSRI(),numMonte=1000)
sim2 = test_convergence(dts,prob,SOSRI2(),numMonte=1000)

p2 = plot(dts,[sim1.errors[:l2],sim2.errors[:l2]],
          title="Linear Convergence Tests",label=["SOSRI" "SOSRI2"],
          ylabel = "Error", xlabel="dt",
          xtickfont = font(16, "LM Roman"),titlefont = font(20, "LM Roman"),
          ytickfont = font(16, "LM Roman"),guidefont = font(18, "LM Roman"),
          legendfont = font(14, "LM Roman"),
          yscale=:log10,xscale=:log10,lw=3)

dts = 1./2.^(7:-1:4) #14->7 good plot
ref = dts.^1.5 * 5e-1
plot!(p2,dts,ref,linestyle=:dash,lw=3,label="1.5 Order Reference")

α = 0.1
β = 0.5
ff1 = (t,u) -> β./sqrt.(1+t)
ff2 = (t,u) -> - u./(2*(1+t))
σ2 = (t,u) -> α*β./sqrt.(1+t)
prob = SplitSDEProblem(ff1,ff2,σ2,1.,(0.0,1.0))
(p::typeof(prob.f))(::Type{Val{:analytic}},t,u0,W) = u0./sqrt.(1+t) + β*(t+α*W)./sqrt.(1+t)

dts = 1./2.^(10:-1:2) #14->7 good plot
sim = test_convergence(dts,prob,RackKenCarp(),numMonte=Int(1e3))

p3 = plot(dts,sim.errors[:l2],
          title="IMEX Convergence Tests",label="IMEX SKenCarp",
          ylabel = "Error", xlabel="dt",
          xtickfont = font(16, "LM Roman"),titlefont = font(20, "LM Roman"),
          ytickfont = font(16, "LM Roman"),guidefont = font(18, "LM Roman"),
          legendfont = font(14, "LM Roman"),
          yscale=:log10,xscale=:log10,lw=3)

ref = dts.^2 * 1e-4
plot!(p3,dts,ref,linestyle=:dash,lw=3,label="2nd Order Reference")

plot(p1,p2,p3,layout=(3,1),size=(600,800))

savefig("convergence.png")
savefig("convergence.pdf")
