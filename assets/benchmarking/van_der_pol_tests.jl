using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
using ParameterizedFunctions, DiffEqMonteCarlo, Base.Test, OrdinaryDiffEq
using Plots; gr()
using BenchmarkTools

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
p0 = plot(sol0,plotdensity=20000,denseplot=true,ylims=[-10,10],title="ODE Solution")

@time sol1 =solve(prob,SOSRA();abstol=1,reltol=1/2^1)
p1 = plot(sol1,plotdensity=20000,denseplot=true,ylims=[-10,10],title="High tolerance SOSRA Solution")
@time sol2 =solve(prob,SOSRA2();abstol=1/2,reltol=1/2^2)
plot(sol2,plotdensity=20000,denseplot=true,ylims=[-10,10])
@time sol3 =solve(prob,SRA1();abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(prob,SRA3();abstol=1/2^3,reltol=1/2^3)
p2 = plot(sol4,plotdensity=20000,denseplot=true,ylims=[-10,10],title="High Tolerance SRA3 Solution")
@time sol42 =solve(prob,SRA3();abstol=1/2^6,reltol=1/2^4)
p3 = plot(sol42,plotdensity=20000,denseplot=true,ylims=[-10,10],title="Low Tolerance SRA3 Solution")
offset = 0

plot(p0,p1,p2,p3,layout=grid(2,2))

######

#### Setup
@time sol1 =solve(prob,SOSRA();abstol=10,reltol=1/2^1)
@time sol2 =solve(prob,SOSRA2();abstol=1/2,reltol=1/2^2)
@time sol3 =solve(prob,SRA1();abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(prob,SRA3();abstol=1/2^6,reltol=1/2^3)
@time sol5 =solve(prob,RackKenCarp();abstol=1,reltol=1)

#### Monte Carlo
@time sol1 =solve(monte_prob,SOSRA();num_monte=N,abstol=10,reltol=1/2^1,verbose=false,save_everystep=false)
@time sol2 =solve(monte_prob,SOSRA2();num_monte=N,abstol=1,reltol=1/2^2,verbose=false,save_everystep=false)
@time sol3 =solve(monte_prob,SRA1();num_monte=N,abstol=1/2^6,reltol=1/2^(5),verbose=false,save_everystep=false)
@time sol4 =solve(monte_prob,SRA3();num_monte=N,abstol=1/2^6,reltol=1/2^3,verbose=false,save_everystep=false)
@time sol5 =solve(monte_prob,RackKenCarp();num_monte=N,abstol=10,reltol=1/2,verbose=false,save_everystep=false)

for i in 1:5
  sol1 =solve(monte_prob,SRA1();num_monte=N,abstol=10,reltol=1/2^(2),verbose=false)
  sol2 =solve(monte_prob,SRA3();num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  sol3 =solve(monte_prob,SOSRA();num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  sol4 =solve(monte_prob,SOSRA2();num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  fails[i,1] = N - sum([sol1[i].retcode == :Success for i in 1:N])
  fails[i,2] = N - sum([sol2[i].retcode == :Success for i in 1:N])
  fails[i,3] = N - sum([sol3[i].retcode == :Success for i in 1:N])
  fails[i,4] = N - sum([sol4[i].retcode == :Success for i in 1:N])
  times[i,1] = sol1.elapsedTime
  times[i,2] = sol2.elapsedTime
  times[i,3] = sol3.elapsedTime
  times[i,4] = sol4.elapsedTime
end

@benchmark sol =solve(prob,SOSRA();abstol=10,reltol=1/2^1,timeseries_steps=10)
@benchmark sol =solve(prob,SOSRA2();abstol=1/2^1,reltol=1/2^2,timeseries_steps=10)
@benchmark sol2 =solve(prob,SRA1();abstol=1/2^6,reltol=1/2^5,timeseries_steps=10)
@benchmark sol2 =solve(prob,SRA3();abstol=1/2^3,reltol=1/2^3,timeseries_steps=10)

#######

# Eigenvalue estimates (put inside of the perform_step)

integrator.f(@muladd(t + c₀[2]*dt),H0[2],ftmp)
integrator.f(@muladd(t + c₀[3]*dt),H0[3],gtmp)
eig_est = dt*(norm(ftmp-gtmp)/norm(H0[2]-H0[3]))/5
eig_est > 5 && println("$t")
