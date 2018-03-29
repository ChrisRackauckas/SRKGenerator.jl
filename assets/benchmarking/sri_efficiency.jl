using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
using ParameterizedFunctions, DiffEqMonteCarlo, Base.Test, OrdinaryDiffEq
using Plots; pyplot()
using BenchmarkTools
names = ["EM","RKMil","SRIW1","SRIW2","SOSRI","SOSRI2"]

################################################################################

prob = prob_sde_2Dlinear

reltols = 1.0./10.0.^(1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) + 2))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) + 2))
          Dict(:alg=>SRIW1())
          Dict(:alg=>SRIW2())
          Dict(:alg=>SOSRI())
          Dict(:alg=>SOSRI2())
          ]

wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns=1000,names=names,error_estimate=:l2)

p1 = plot(wp,legend=false,title="Linear Strong Error",
          xtickfont = font(16, "LM Roman"),titlefont = font(20, "LM Roman"),
          ytickfont = font(16, "LM Roman"),guidefont = font(18, "LM Roman"))

wp2= WorkPrecisionSet(prob,abstols,reltols,setups;numruns=1000,names=names,error_estimate=:weak_final)

sample_size = Int[10;1e2;1e3;1e4;1e5]
test_dt = 1e-2
appxsol_setup = Dict(:alg=>SOSRA(),:abstol=>1e-4,:reltol=>1e-4)
se1 = get_sample_errors(prob,test_dt,appxsol_setup=appxsol_setup,
          parallel_type = :threads, numruns=sample_size, std_estimation_runs = Int(1e3))

p2 = plot(wp2;legend=false,title="Linear Weak Error",
          xtickfont = font(16, "LM Roman"),titlefont = font(20, "LM Roman"),
          ytickfont = font(16, "LM Roman"),guidefont = font(18, "LM Roman"),
          plot_sample_error = false)
times = [wp2[i].times for i in 1:length(wp)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!(p2,[se1[2];se1[2]],times,color=:red,linestyle=:dash,label="Sample Error: 100",lw=3)
plot!(p2,[se1[end];se1[end]],times,color=:orange,linestyle=:dash,label="Sample Error: 10000",lw=3)

### Lotka

f = @ode_def LotkaVolterraTest begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=1.0 c=3.0 d=1.0

function g(t,u,du)
  @. du = 0.01u
end
u0 = [1.0;1.0]
tspan = (0.0,10.0)
prob = SDEProblem(f,g,u0,tspan);

reltols = 1.0./4.0.^(2:4)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [Dict(:alg=>SRIW1())
          Dict(:alg=>EM(),:dts=>1.0./12.0.^((1:length(reltols)) + 1.5))
          Dict(:alg=>RKMil(),:dts=>1.0./12.0.^((1:length(reltols)) + 1.5))
          Dict(:alg=>SRIW2())
          Dict(:alg=>SOSRI())
          Dict(:alg=>SOSRI2())
          ]
test_dt = 1/10^2
appxsol_setup = Dict(:alg=>SRIW1(),:abstol=>1e-4,:reltol=>1e-4)
wp3 = WorkPrecisionSet(prob,abstols,reltols,setups,test_dt;
                                     names = names,
                                     verbose=false,save_everystep=false,
                                     parallel_type = :threads,
                                     appxsol_setup = appxsol_setup,
                                     numruns_error=100,error_estimate=:final)
plot(wp3)
p3 = plot(wp3,legend=false,title="Lotka-Volterra Strong Error",
          xtickfont = font(16, "LM Roman"),titlefont = font(20, "LM Roman"),
          ytickfont = font(16, "LM Roman"),guidefont = font(18, "LM Roman"))

wp4 = WorkPrecisionSet(prob,abstols,reltols,setups,test_dt;
                                     names = names,
                                     verbose=false,save_everystep=false,
                                     parallel_type = :threads,
                                     appxsol_setup = appxsol_setup,
                                     numruns_error=100,error_estimate=:weak_final)
plot(wp4)
sample_size = Int[10;1e2;1e3;1e4;1e5]
test_dt = 1e-2
appxsol_setup = Dict(:alg=>SOSRA(),:abstol=>1e-4,:reltol=>1e-4)
se2 = get_sample_errors(prob,test_dt,appxsol_setup=appxsol_setup,
          parallel_type = :threads, numruns=sample_size, std_estimation_runs = Int(1e3))

p4 = plot(wp4;legend=:topleft,title="Lotka-Volterra Weak Error",
          plot_sample_error = false,
          xtickfont = font(16, "LM Roman"),titlefont = font(20, "LM Roman"),
          ytickfont = font(16, "LM Roman"),guidefont = font(18, "LM Roman"),
          legendfont = font(14, "LM Roman"))
times = [wp4[i].times for i in 1:length(wp4)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!(p4,[se2[2];se2[2]],times,color=:red,linestyle=:dash,label="Sample Error: 100",lw=3)
plot!(p4,[se2[end];se2[end]],times,color=:orange,linestyle=:dash,label="Sample Error: 10000",lw=3)

################################################################################

### Final Plot

plot(p1,p2,p3,p4,layout=(2,2),size=(1200,800))

savefig("SRI_efficiency.png")
savefig("SRI_efficiency.pdf")
