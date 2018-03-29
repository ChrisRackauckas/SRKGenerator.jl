using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
using ParameterizedFunctions, DiffEqMonteCarlo, Base.Test, OrdinaryDiffEq
using Plots; pyplot()
using BenchmarkTools, DiffEqNoiseProcess

################################################################################

### Van Der Pol

van = @ode_def VanDerPol2 begin
  dy = μ*((1-x^2)*y - x)
  dx = 1*y
end μ

σ = (du,u,p,t) -> begin
  for i in 1:2
    du[i] = 3.0 #Additive
  end
end

ode_prob = ODEProblem(van,[0;2.],(0.0,6.3),(1e5,))
prob = SDEProblem(van,σ,[0;2.],(0.0,6.3),(1e5,))
monte_prob = MonteCarloProblem(prob)
N = 100; fails = zeros(5,4); times = zeros(5,4)
offset = 0

#####

@time sol0 =solve(ode_prob,Tsit5();abstol=1e-6,reltol=1e-3)
p0 = plot(sol0,plotdensity=20000,denseplot=true,ylims=[-10,10],title="ODE Solution",
          xtickfont = font(16, "Latin Modern Roman"),
          titlefont = font(20, "Latin Modern Roman"),
          ytickfont = font(16, "Latin Modern Roman"),
          guidefont = font(18, "Latin Modern Roman"),
          legendfont = font(14, "Latin Modern Roman"))

srand(101)
test_dt = 0.0001
t = prob.tspan[1]:test_dt:prob.tspan[2]
brownian_values = cumsum([[zeros(size(prob.u0))];[sqrt(test_dt)*randn(size(prob.u0)) for i in 1:length(t)-1]])
brownian_values2 = cumsum([[zeros(size(prob.u0))];[sqrt(test_dt)*randn(size(prob.u0)) for i in 1:length(t)-1]])
np = NoiseGrid(t,brownian_values,brownian_values2)
fixed_noise_prob = SDEProblem(van,σ,[0;2.],(0.0,6.3),(1e5,),noise=np)

@time sol1 =solve(fixed_noise_prob,SOSRA();abstol=1,reltol=1/2^1)
p1 = plot(sol1,plotdensity=20000,denseplot=true,ylims=[-10,10],
          title="High tolerance SOSRA",
          legend = false,
          xtickfont = font(16, "Latin Modern Roman"),
          titlefont = font(20, "Latin Modern Roman"),
          ytickfont = font(16, "Latin Modern Roman"),
          guidefont = font(18, "Latin Modern Roman"),
          legendfont = font(14, "Latin Modern Roman"))

#@time sol1 =solve(prob,ImplicitEM();dt=1/1000000) # Errors

@time sol2 =solve(fixed_noise_prob,SOSRA2();abstol=1/2,reltol=1/2^2)
plot(sol2,plotdensity=20000,denseplot=true,ylims=[-10,10])
#@time sol3 =solve(prob,SRA1();abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(fixed_noise_prob,SRA3();abstol=1/2^3,reltol=1/2^3)
p2 = plot(sol4,plotdensity=20000,denseplot=true,ylims=[-10,10],
          title="High Tolerance SRA3",
          legend = false,
          xtickfont = font(16, "Latin Modern Roman"),
          titlefont = font(20, "Latin Modern Roman"),
          ytickfont = font(16, "Latin Modern Roman"),
          guidefont = font(18, "Latin Modern Roman"),
          legendfont = font(14, "Latin Modern Roman"))
@time sol42 =solve(fixed_noise_prob,SRA3();abstol=1/2^6,reltol=1/2^4)
p3 = plot(sol42,plotdensity=20000,denseplot=true,ylims=[-10,10],
          title="Low Tolerance SRA3",
          legend = false,
          xtickfont = font(16, "Latin Modern Roman"),
          titlefont = font(20, "Latin Modern Roman"),
          ytickfont = font(16, "Latin Modern Roman"),
          guidefont = font(18, "Latin Modern Roman"),
          legendfont = font(14, "Latin Modern Roman"))

@time sol52 =solve(fixed_noise_prob,SKenCarp();abstol=16,reltol=16)
p5 = plot(sol52,plotdensity=20000,denseplot=true,ylims=[-10,10])

diff(sol52.t)

plot(p0,p1,p2,p3,layout=grid(2,2),size=(1200,800))

savefig("additive_van_der_pol.png")
savefig("additive_van_der_pol.pdf")

######

#### Setup
@time sol1 =solve(prob,SOSRA();abstol=10,reltol=1/2^1)
@time sol2 =solve(prob,SOSRA2();abstol=1/2,reltol=1/2^2)
@time sol3 =solve(prob,SRA1();abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(prob,SRA3();abstol=1/2^6,reltol=1/2^3)
@time sol5 =solve(prob,SKenCarp();abstol=1,reltol=1)

#### Estimates

prob_quick = SDEProblem(VanDerPol2(μ=1e5),σ,[0;2.],(0.0,1.0))
@time sol6 =solve(prob_quick,EM();dt=5e-7) # Unstable
@time sol6 =solve(prob_quick,EM();dt=1e-7) # Stable?

42.087122 *6.3 * 100 # 7 hours 21 minutes....

@time sol6 =solve(prob,EM();dt=1e-7) # Stable?

346.117152 * 100 # 9 hours 36 minutes 52 seconds

@time sol6 =solve(prob,EM();dt=5e-8) # Stable?

@time sol6 =solve(prob_quick,ImplicitEM();dt=1e-7) # Unstable
@time sol6 =solve(prob_quick,ImplicitEM();dt=5e-8) # Stable?

363.668494 * 6.3 * 100 # 2 days 15 hours 38 minutes and 31 seconds...

#### Monte Carlo
println("SOSRA time")
@time sol1 =solve(monte_prob,SOSRA();num_monte=N,abstol=10,reltol=1/2^1,verbose=false,save_everystep=false)
println("SOSRA2 time")
@time sol2 =solve(monte_prob,SOSRA2();num_monte=N,abstol=1,reltol=1/2^2,verbose=false,save_everystep=false)
println("SRA1 time")
@time sol3 =solve(monte_prob,SRA1();num_monte=N,abstol=1/2^6,reltol=1/2^(5),verbose=false,save_everystep=false)
println("SRA3 time")
@time sol4 =solve(monte_prob,SRA3();num_monte=N,abstol=1/2^6,reltol=1/2^3,verbose=false,save_everystep=false)
println("SKenCarp time")
@time sol5 =solve(monte_prob,SKenCarp();num_monte=N,abstol=16,reltol=16,verbose=false,save_everystep=false)

sum(sol.retcode == :Success for sol in sol5)

@time sol5 =solve(monte_prob,SKenCarp();num_monte=N,abstol=4,reltol=4,verbose=false,save_everystep=false)

sum(sol.retcode == :Success for sol in sol5)

@time sol5 =solve(monte_prob,SKenCarp(extrapolant=:trivial);num_monte=N,
                  abstol=4,reltol=4,verbose=false,save_everystep=false)

sum(sol.retcode == :Success for sol in sol5)

println("EM time")
# dt = 1e-7 has 15 failures
@time sol5 =solve(monte_prob,EM();num_monte=N,dt=1e-8,save_everystep=false)

sum(sol.retcode == :Success for sol in sol5)

#######

# Eigenvalue estimates (put inside of the perform_step)

integrator.f(@muladd(t + c₀[2]*dt),H0[2],ftmp)
integrator.f(@muladd(t + c₀[3]*dt),H0[3],gtmp)
eig_est = dt*(norm(ftmp-gtmp)/norm(H0[2]-H0[3]))/5
eig_est > 5 && println("$t")
