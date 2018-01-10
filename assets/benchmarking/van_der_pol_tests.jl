using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
using ParameterizedFunctions, DiffEqMonteCarlo, Base.Test, OrdinaryDiffEq
using Plots; pyplot()
using BenchmarkTools, DiffEqNoiseProcess

################################################################################

### Van Der Pol

van = @ode_def VanDerPol2 begin
  dy = μ*((1-x^2)*y - x)
  dx = 1*y
end μ=>1.

σ = (t,u,du) -> begin
  for i in 1:2
    du[i] = 3.0 #Additive
  end
end

ode_prob = ODEProblem(VanDerPol2(μ=1e5),[0;2.],(0.0,6.3))
prob = SDEProblem(VanDerPol2(μ=1e5),σ,[0;2.],(0.0,6.3))
monte_prob = MonteCarloProblem(prob)
N = 100; fails = zeros(5,4); times = zeros(5,4)
offset = 0

#####

@time sol0 =solve(ode_prob,Tsit5();abstol=1e-6,reltol=1e-3)
p0 = plot(sol0,plotdensity=20000,denseplot=true,ylims=[-10,10],title="ODE Solution",
          xtickfont = font(16, "Ariel"),titlefont = font(20, "Ariel"),
          ytickfont = font(16, "Ariel"),guidefont = font(18, "Ariel"),
          legendfont = font(14, "Ariel"))

srand(101)
test_dt = 0.0001
t = prob.tspan[1]:test_dt:prob.tspan[2]
brownian_values = cumsum([[zeros(size(prob.u0))];[sqrt(test_dt)*randn(size(prob.u0)) for i in 1:length(t)-1]])
brownian_values2 = cumsum([[zeros(size(prob.u0))];[sqrt(test_dt)*randn(size(prob.u0)) for i in 1:length(t)-1]])
np = NoiseGrid(t,brownian_values,brownian_values2)
fixed_noise_prob = SDEProblem(VanDerPol2(μ=1e5),σ,[0;2.],(0.0,6.3),noise=np)

@time sol1 =solve(fixed_noise_prob,SOSRA();abstol=1,reltol=1/2^1,seed=5)
p1 = plot(sol1,plotdensity=20000,denseplot=true,ylims=[-10,10],title="High tolerance SOSRA",
          legend = false,
          xtickfont = font(16, "Ariel"),titlefont = font(20, "Ariel"),
          ytickfont = font(16, "Ariel"),guidefont = font(18, "Ariel"),
          legendfont = font(14, "Ariel"))

#@time sol1 =solve(prob,ImplicitEM();dt=1/1000000) # Errors

@time sol2 =solve(fixed_noise_prob,SOSRA2();abstol=1/2,reltol=1/2^2)
plot(sol2,plotdensity=20000,denseplot=true,ylims=[-10,10])
#@time sol3 =solve(prob,SRA1();abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(fixed_noise_prob,SRA3();abstol=1/2^3,reltol=1/2^3)
p2 = plot(sol4,plotdensity=20000,denseplot=true,ylims=[-10,10],title="High Tolerance SRA3",
          legend = false,
          xtickfont = font(16, "Ariel"),titlefont = font(20, "Ariel"),
          ytickfont = font(16, "Ariel"),guidefont = font(18, "Ariel"),
          legendfont = font(14, "Ariel"))
@time sol42 =solve(fixed_noise_prob,SRA3();abstol=1/2^6,reltol=1/2^4)
p3 = plot(sol42,plotdensity=20000,denseplot=true,ylims=[-10,10],title="Low Tolerance SRA3",
          legend = false,
          xtickfont = font(16, "Ariel"),titlefont = font(20, "Ariel"),
          ytickfont = font(16, "Ariel"),guidefont = font(18, "Ariel"),
          legendfont = font(14, "Ariel"))

#@time sol52 =solve(fixed_noise_prob,RackKenCarp();abstol=1/2^6,reltol=1/2^6)
#p5 = plot(sol52,plotdensity=20000,denseplot=true,ylims=[-10,10])

plot(p0,p1,p2,p3,layout=grid(2,2),size=(1200,800))

savefig("additive_van_der_pol.png")
savefig("additive_van_der_pol.pdf")

######

#### Setup
@time sol1 =solve(prob,SOSRA();abstol=10,reltol=1/2^1)
@time sol2 =solve(prob,SOSRA2();abstol=1/2,reltol=1/2^2)
@time sol3 =solve(prob,SRA1();abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(prob,SRA3();abstol=1/2^6,reltol=1/2^3)
@time sol5 =solve(prob,RackKenCarp();abstol=1,reltol=1)

#### Monte Carlo
println("SOSRA time")
@time sol1 =solve(monte_prob,SOSRA();num_monte=N,abstol=10,reltol=1/2^1,verbose=false,save_everystep=false)
println("SOSRA2 time")
@time sol2 =solve(monte_prob,SOSRA2();num_monte=N,abstol=1,reltol=1/2^2,verbose=false,save_everystep=false)
println("SRA1 time")
@time sol3 =solve(monte_prob,SRA1();num_monte=N,abstol=1/2^6,reltol=1/2^(5),verbose=false,save_everystep=false)
println("SRA3 time")
@time sol4 =solve(monte_prob,SRA3();num_monte=N,abstol=1/2^6,reltol=1/2^3,verbose=false,save_everystep=false)
println("RackKenCarp time")
@time sol5 =solve(monte_prob,RackKenCarp();num_monte=N,abstol=1/2^6,reltol=1/2^6,verbose=false,save_everystep=false)

#######

# Eigenvalue estimates (put inside of the perform_step)

integrator.f(@muladd(t + c₀[2]*dt),H0[2],ftmp)
integrator.f(@muladd(t + c₀[3]*dt),H0[3],gtmp)
eig_est = dt*(norm(ftmp-gtmp)/norm(H0[2]-H0[3]))/5
eig_est > 5 && println("$t")
