using FFTW

using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
gr()

prob = prob_sde_additive
prob.tspan = (0.0,1.0)

reltols = 1.0./10.0.^(-2:4)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) - 2))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) - 2))
          Dict(:alg=>SRA(tableau=StochasticDiffEq.constructSRA1()))
          Dict(:alg=>SRA(tableau=StochasticDiffEq.constructSRA2()))
          Dict(:alg=>SRA(tableau=StochasticDiffEq.constructSRA3()))
          Dict(:alg=>SRA(tableau=StochasticDiffEq.constructSOSRA()))
          Dict(:alg=>SRA(tableau=StochasticDiffEq.constructSOSRA2()))
          ]
names = ["EM","RKMil","SRA1","SRA2","SRA3","SOSRA","SOSRA2"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns=10000,names=names,error_estimate=:l2)

dts = 1./2.^(10:-1:2) #14->7 good plot
using Base.Test
sim4 = test_convergence(dts,prob,SRA(tableau=StochasticDiffEq.constructSOSRA()),numMonte=Int(1e1))
@test abs(sim4.𝒪est[:final]-2) < 0.3
sim5 = test_convergence(dts,prob,SRA(tableau=StochasticDiffEq.constructSOSRA2()),numMonte=Int(1e1))
@test abs(sim5.𝒪est[:final]-2) < 0.3

using Plots; pyplot()
p1 = plot(wp,title="Work-Precision Plot")
p2 = plot(sim4,title="Convergence Tests",label=["SOSRA l∞" "SOSRA weak final error" "SOSRA final error" "SOSRA l2"])
plot!(p2,sim5,label=["SOSRA2 l∞" "SOSRA2 weak final" "SOSRA2 final" "SOSRA2 l2"])
plot(p1,p2,layout=(2,1))


using FFTW
using StochasticDiffEq, DiffEqProblemLibrary, ParameterizedFunctions
srand(100)

#=
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
=#

using DiffEqMonteCarlo
monte_prob = MonteCarloProblem(prob_sde_lorenz)
N = 100; fails = zeros(5,4); times = zeros(5,4)
offset = 0
for i in 1:5
  sol1 =solve(monte_prob,SRA1();num_monte=N,abstol=10,reltol=1/2^(2),verbose=false)
  sol2 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSRA3());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  sol3 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSOSRA());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  sol4 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  fails[i,1] = N - sum([sol1[i].retcode == :Success for i in 1:N])
  fails[i,2] = N - sum([sol2[i].retcode == :Success for i in 1:N])
  fails[i,3] = N - sum([sol3[i].retcode == :Success for i in 1:N])
  fails[i,4] = N - sum([sol4[i].retcode == :Success for i in 1:N])
  times[i,1] = sol1.elapsedTime
  times[i,2] = sol2.elapsedTime
  times[i,3] = sol3.elapsedTime
  times[i,4] = sol4.elapsedTime
end

van = @ode_def_noinvhes VanDerPol2 begin
  dy = μ*((1-x^2)*y - x)
  dx = 1*y
end μ=>1.

σ = (t,u,du) -> begin
  for i in 1:2
    du[i] = 3.0 #Additive
  end
end

using OrdinaryDiffEq
prob = ODEProblem(VanDerPol2(μ=1e5),[0;2.],(0.0,6.3))
@time sol0 =solve(prob,Tsit5();abstol=1e-6,reltol=1e-3)
p0 = plot(sol0,plotdensity=20000,denseplot=true,ylims=[-10,10])

prob = SDEProblem(VanDerPol2(μ=1e5),σ,[0;2.],(0.0,6.3))
monte_prob = MonteCarloProblem(prob)
N = 100; fails = zeros(5,4); times = zeros(5,4)
offset = 0

@time sol1 =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA());abstol=1,reltol=1/2^1)
p1 = plot(sol1,plotdensity=20000,denseplot=true,ylims=[-10,10])
@time sol2 =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());abstol=1/2,reltol=1/2^2)
plot(sol2,plotdensity=20000,denseplot=true,ylims=[-10,10])
@time sol3 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA1());abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA3());abstol=1/2^3,reltol=1/2^3)
p2 = plot(sol4,plotdensity=20000,denseplot=true,ylims=[-10,10])
@time sol42 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA3());abstol=1/2^6,reltol=1/2^4)
p3 = plot(sol42,plotdensity=20000,denseplot=true,ylims=[-10,10])
offset = 0

plot(p0,p1,p2,p3,layout=grid(2,2))

N = 100

#### Setup
@time sol1 =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA());abstol=10,reltol=1/2^1)
@time sol2 =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());abstol=1/2,reltol=1/2^2)
@time sol3 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA1());abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA3());abstol=1/2^6,reltol=1/2^3)
#### Monte Carlo
@time sol1 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSOSRA());num_monte=N,abstol=10,reltol=1/2^1,verbose=false,save_everystep=false)
@time sol2 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());num_monte=N,abstol=1,reltol=1/2^2,verbose=false,save_everystep=false)
@time sol3 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSRA1());num_monte=N,abstol=1/2^6,reltol=1/2^(5),verbose=false,save_everystep=false)
@time sol4 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSRA3());num_monte=N,abstol=1/2^6,reltol=1/2^3,verbose=false,save_everystep=false)



for i in 1:5
  sol1 =solve(monte_prob,SRA1();num_monte=N,abstol=10,reltol=1/2^(2),verbose=false)
  sol2 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSRA3());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  sol3 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSOSRA());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  sol4 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  fails[i,1] = N - sum([sol1[i].retcode == :Success for i in 1:N])
  fails[i,2] = N - sum([sol2[i].retcode == :Success for i in 1:N])
  fails[i,3] = N - sum([sol3[i].retcode == :Success for i in 1:N])
  fails[i,4] = N - sum([sol4[i].retcode == :Success for i in 1:N])
  times[i,1] = sol1.elapsedTime
  times[i,2] = sol2.elapsedTime
  times[i,3] = sol3.elapsedTime
  times[i,4] = sol4.elapsedTime
end
using BenchmarkTools
@benchmark sol =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA());abstol=10,reltol=1/2^1,timeseries_steps=10)
@benchmark sol =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());abstol=1/2^1,reltol=1/2^2,timeseries_steps=10)
@benchmark sol2 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA1());abstol=1/2^6,reltol=1/2^5,timeseries_steps=10)
@benchmark sol2 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA3());abstol=1/2^3,reltol=1/2^3,timeseries_steps=10)





srand(220)
prob2 = SDEProblem(VanDerPol2(μ=1e5),σ,[0;2.],(0.0,6.3))
@time sol2 = solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());abstol=1/2,reltol=1/2^2)
p1 = plot(sol2,plotdensity=20000,denseplot=true,ylims=[-10,10],title="Trajectory 1")
ts1 =[0.9764368351990603,0.9775967457127613,0.978859704361288,0.9798533139430902,0.9817240075118023,0.9825000582530085,0.9836091918077459,0.9852152561862988,0.9860667766563265,0.9865652370522721,1.1878855302678386,1.1903249867063068,1.1904526365864387,1.1934139873164036,1.1934740025075103,1.1941047010314014,1.640902086834412,1.6450022275786294,2.0570569909698406,2.0640351051104204,2.0676863230571736,2.4042117334572373,2.405057062025887,2.463688542881055,2.467842908983563,2.469948048177236,2.4700286399169595,2.470221304544736,2.47059392806094,2.512628958932553,2.5329847996389763,2.5335746404484643,2.6427664310243655,2.649598618156501,2.654805240550594,2.6550483975565147,2.655456731954062,2.8310123878389524,2.831209926303201,2.831922838390023,2.851968951521406,2.8534434995619424,2.8538002183021245,2.871439832658957,2.871984204663848,2.873047992944743,2.8802745447585423,2.880451323539189,2.881657926850309,2.8830851867846228,2.883225979216535,2.8848311965671125,2.88498058657672,2.8851995146372666,2.8856303310061504,2.8863329966360056,3.1605943805112857,3.1620751699111143,3.1625284900532526,3.16436315156008,6.273517199217812,6.27653837434436,6.276711322636875,6.2770918582267035,6.279750375845481,6.27979651640385,6.2799655792381905,6.280427537384802,6.280744810844847]
scatter!(p1,ts1,markersize=10,zeros(length(ts)),ylims=[-10,10],label="Stiffness Detection Positive")

srand(230)
prob2 = SDEProblem(VanDerPol2(μ=1e5),σ,[0;2.],(0.0,6.3))
@time sol22 =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());abstol=1/2,reltol=1/2^2)
p2 = plot(sol22,plotdensity=20000,denseplot=true,ylims=[-10,10],legend=false,title="Trajectory 2")
ts2 = [0.046616479937447504,0.04745354112633021,0.06871511983609467,0.06974319280424324,0.06979634692312843,0.07015975175606998,0.07043055696092873,0.07149372457591815,0.07184149050047958,0.11595416067865726,0.11615787740732275,0.11850614200604724,0.11886154580024354,0.13757344004569674,0.13764221324043185,0.13771958308450888,0.14203677115919802,0.14266447662479878,0.14292081554544428,0.14323620011110305,0.1435377758333027,0.14425595378611122,0.14434224005381127,0.17696091173191408,0.18025473413506166,0.18140247158672124,0.18195479903110687,0.18201724449653697,0.4398821665826248,0.44014030978420693,0.44199588257597716,0.442232743993479,0.44246149150761477,0.44253103664119264,0.44307943601925154,0.44317203608265054,0.4491050689774769,0.4492283721485624,0.44941505062693865,0.44948890702568434,0.45005860645944545,0.45231144374012994,0.45761542564573865,0.45767481697207246,0.4578401894124175,0.45943736712601163,0.459665866299597,0.4599729359551946,0.46011411068843383,0.46190153920906185,0.46198949957755453,0.46216457957437784,0.4623960909833458,0.46297288186350183,0.4633002984072603,0.46363957443920706,0.6925709260939786,0.7057327295805634,0.7090100792436287,0.7132090910891157,0.7219379557383332,0.7223144560985898,0.8714417944277604,0.8725889128823758,0.8737946688299157,0.8742268192449316,0.9142412113410381,0.9152813070192135,0.9154827360976983,0.915767540244774,0.9168835057011164,0.9171294998097814,0.9172852194322852,1.385361635489688,1.3869127528869964,1.3873005444593085,1.3882147961149047,1.3884991283537416,1.3886661902952382,1.5276285785703287,1.5279104366224654,1.5833909768748715,1.6660901154670178,1.6672269291342856,1.667389680031175,1.669550006357094,1.6708890496710056,2.002353939313724,2.002683993599299,2.031310640608775,5.089757984832253,5.091822398881064,5.092750842009686,5.103639091088352,5.10369159102625,5.10380009477738,5.151588708876289,5.152010247921434,5.152678807586772,5.189220144464494,5.189390131639592,5.189469928767664,5.189555827624455,5.7516101240775805,5.752134323342594,5.7692714454091405,5.7694688725258905,5.774408096395464,5.777630799644512,5.783874527607929,5.7981920981101265,5.87508098463563,5.878051348837662,5.878323021334675,5.879359699316764,5.911856704666698,5.912553111298001,5.914149208822204,5.916747938470164,5.917008076816422,5.917094626902111,5.917740025854637,5.920254856715425,5.921289658564189,5.924406920990269,5.927078735345617,5.9271893916654905,5.927614385097295,5.9285028751815405,5.929360544531681,5.929707355108277,5.932918067982112,5.933048187098134,5.938036986603212,5.948059590991143,5.979103304328124,6.197397040129443,6.241729179166541,6.242067055580289,6.242124698040826,6.244729810731595]
scatter!(p2,ts2,markersize=10,zeros(length(ts)),ylims=[-10,10],label="Stiffness Detection Positive")

plot(p1,p2,layout=grid(2,1))

for i in 1:5
  sol1 =solve(monte_prob,SRA1();num_monte=N,abstol=10,reltol=1/2^(2),verbose=false)
  sol2 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSRA3());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  sol3 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSOSRA());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  sol4 =solve(monte_prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());num_monte=N,abstol=10,reltol=1/2^(offset+i),verbose=false)
  fails[i,1] = N - sum([sol1[i].retcode == :Success for i in 1:N])
  fails[i,2] = N - sum([sol2[i].retcode == :Success for i in 1:N])
  fails[i,3] = N - sum([sol3[i].retcode == :Success for i in 1:N])
  fails[i,4] = N - sum([sol4[i].retcode == :Success for i in 1:N])
  times[i,1] = sol1.elapsedTime
  times[i,2] = sol2.elapsedTime
  times[i,3] = sol3.elapsedTime
  times[i,4] = sol4.elapsedTime
end
using BenchmarkTools
@benchmark sol =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA());abstol=10,reltol=1/2^1,timeseries_steps=10)
@benchmark sol =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());abstol=1/2^1,reltol=1/2^2,timeseries_steps=10)
@benchmark sol2 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA1());abstol=1/2^6,reltol=1/2^5,timeseries_steps=10)
@benchmark sol2 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA3());abstol=1/2^3,reltol=1/2^3,timeseries_steps=10)


using OrdinaryDiffEq
prob = ODEProblem(VanDerPol2(μ=1e5),[0;2.],(0.0,6.3))
@time sol0 =solve(prob,Tsit5();abstol=1e-6,reltol=1e-3)
p0 = plot(sol0,plotdensity=20000,denseplot=true,ylims=[-10,10],title="ODE Solution")

prob = SDEProblem(VanDerPol2(μ=1e5),σ,[0;2.],(0.0,6.3))
monte_prob = MonteCarloProblem(prob)
N = 100; fails = zeros(5,4); times = zeros(5,4)
offset = 0

@time sol1 =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA());abstol=1,reltol=1/2^1)
p1 = plot(sol1,plotdensity=20000,denseplot=true,ylims=[-10,10],title="High tolerance SOSRA Solution")
@time sol2 =solve(prob,SRA(tableau=StochasticDiffEq.constructSOSRA2());abstol=1/2,reltol=1/2^2)
plot(sol2,plotdensity=20000,denseplot=true,ylims=[-10,10])
ts = [0.1287417376231843,0.12885855628048268,0.22313330634159534,0.22557815771064949,0.22607078639208922,0.4270211579553978,0.42766124199177236,0.42824255808757333,0.4284754913204088,0.5252536670748583,1.7430487638156489,1.7465124731868622,1.7467417097520845,1.8342414009259005,2.8532223352688746,3.005847515916775,3.006288440218983,3.007523685449149,3.0077656658549503,3.3100578839942596,3.310148528784323,3.312720190907894,3.3135251461147446,3.3136932079970496,3.4230363009002405,3.5434929094675107,3.5436814003217227,3.5438496007217983,3.5445772536082387,3.6642380174205673,3.664389803952252,3.6656366384033183,3.665834547040817,3.738678399228814,3.740094284628105,3.7407436712099176,3.7578144750777516,3.757941626741623,3.758829321047308,3.900802383866269,3.9015119286233193,5.079780530301039,5.086111739471306,5.087359385591886]
scatter!(ts,zeros(length(ts)),ylims=[-10,10])
@time sol3 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA1());abstol=1/2^7,reltol=1/2^5)
@time sol4 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA3());abstol=1/2^3,reltol=1/2^3)
p2 = plot(sol4,plotdensity=20000,denseplot=true,ylims=[-10,10],title="High Tolerance SRA3 Solution")
@time sol42 =solve(prob,SRA(tableau=StochasticDiffEq.constructSRA3());abstol=1/2^6,reltol=1/2^4)
p3 = plot(sol42,plotdensity=20000,denseplot=true,ylims=[-10,10],title="Low Tolerance SRA3 Solution")
offset = 0

plot(p0,p1,p2,p3,layout=grid(2,2))

integrator.f(@muladd(t + c₀[2]*dt),H0[2],ftmp)
integrator.f(@muladd(t + c₀[3]*dt),H0[3],gtmp)
eig_est = dt*(norm(ftmp-gtmp)/norm(H0[2]-H0[3]))/5
eig_est > 5 && println("$t")


using FFTW
using StochasticDiffEq, Plots, DiffEqDevTools, DiffEqProblemLibrary
pyplot()

prob = prob_sde_2Dlinear
prob.tspan = (0.0,1.0)

reltols = 1.0./10.0.^(1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) + 2))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) + 2))
          Dict(:alg=>SRI())
          Dict(:alg=>SRI(tableau=StochasticDiffEq.constructSRIOpt1()))
          Dict(:alg=>SRI(tableau=StochasticDiffEq.constructSRIOpt2()))
          ]
names = ["EM","RKMil","SRIW1","SOSRI","SOSRI2"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns=1000,names=names,error_estimate=:l2)

using Base.Test
dts = 1./2.^(7:-1:4) #14->7 good plot
srand(100)
prob = prob_sde_2Dlinear
sim1 = test_convergence(dts,prob,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),numMonte=10000)
@test abs(sim1.𝒪est[:final]-1.5) < 0.3
sim2 = test_convergence(dts,prob,SRI(tableau=StochasticDiffEq.constructSRIOpt2()),numMonte=10000)
@test abs(sim2.𝒪est[:final]-1.5) < 0.3

using Plots; pyplot()
p1 = plot(wp,title="Work-Precision Plot")
p2 = plot(sim1,title="Convergence Tests",label=["SOSRI l∞" "SOSRI weak final error" "SOSRI final error" "SOSRI l2"])
plot!(p2,sim2,label=["SOSRI2 l∞" "SOSRI2 weak final error" "SOSRI2 final" "SOSRI2 l2"])
plot(p1,p2,layout=(2,1))
