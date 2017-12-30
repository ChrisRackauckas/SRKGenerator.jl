using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
using ParameterizedFunctions, DiffEqMonteCarlo, Base.Test, OrdinaryDiffEq
using Plots; plotly()
using BenchmarkTools

################################################################################

## Lorenz

f = @ode_def_nohes LorenzSDE begin
  dx = σ*(y-x)
  dy = x*(ρ-z) - y
  dz = x*y - β*z
end σ=>10. ρ=>28. β=>2.66

σ = (t,u,du) -> begin
  for i in 1:3
    du[i] = 10.0 #Additive
  end
end

prob_sde_lorenz = SDEProblem(f,σ,ones(3),(0.0,10.0))

monte_prob = MonteCarloProblem(prob_sde_lorenz)
N = 100; fails = zeros(5,4); times = zeros(5,4)
offset = 0
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
