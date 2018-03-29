using DiffEqBase, OrdinaryDiffEq, StochasticDiffEq, Sundials
using Plots; plotly()

############ Spatial Parameters

const x_start = -100
const x_end   = 400
const y_start = 0
const y_end   = 100
const dx = 5
const dy = 5
const N = Int((x_end-x_start) / dx + 1)
const M = Int((y_end-y_start) / dy +1)
const X = repmat(collect(x_start:dx:x_end)',M,1)
const Y = repmat(collect(y_start:dy:y_end),1,N)

############ RA Parameters

const DRA = 250.46
const β = 1
const kA = 0.002
const kdeg = 50
const kmax = 1
const n = 2
const xf = 400

const kp = 0.1
const gamma = 100
const mon = 3.0
const moff = 0.0013
const rdeg1 = 0.0001
const rdeg2 = 0.0001
const jalpha = 50.0
const jbeta = .10
const bpdeg1 = 0.0001
const d = 1
const e = 1
const Vbp = 1e-3
const Vr = 2e-2
const ϵout = 0.1
const ϵin = 0.1
const ϵR = 0.1

############ Setup Diffusion Matrices
const Ax = Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1])
const Ay = Tridiagonal([1.0 for i in 1:M-1],[-2.0 for i in 1:M],[1.0 for i in 1:M-1])
Ax[2,1] = 2
Ax[end-1,end] = 2
Ay[1,2] = 2
Ay[end,end-1] = 2

Ax ./= dx^2
Ay ./= dy^2

############# Spatial Constants

const Vmax = 20.0
V(x) = ifelse(x>xf-40,Vmax,zero(typeof(Vmax)))
const VRA = V.(X)
no_cyp(x) = ifelse(x<0 || (xf-40<x && x<=xf),true,false)
const NCYP_kmax = no_cyp.(X)

u0 = zeros(M,N,6)

const Ax_cache = zeros(M,N)
const Ay_cache = zeros(M,N)
const diffRA = zeros(M,N)

function ra_gradient(t,u,du)
  RAout = @view u[:,:,1]
  RAin = @view u[:,:,2]
  R = @view u[:,:,3]
  RAR = @view u[:,:,4]
  BP = @view u[:,:,5]
  RABP = @view u[:,:,6]
  dRAout = @view du[:,:,1]
  dRAin = @view du[:,:,2]
  dR = @view du[:,:,3]
  dRAR = @view du[:,:,4]
  dBP = @view du[:,:,5]
  dRABP = @view du[:,:,6]

  A_mul_B!(Ax_cache,RAout,Ax)
  A_mul_B!(Ay_cache,Ay,RAout)

  Threads.@threads for i in eachindex(RAin)
    @inbounds diffRA[i] = DRA*(Ax_cache[i] + Ay_cache[i])
    @inbounds dRAout[i]  = VRA[i] - β*RAout[i] + kp*RAin[i] + diffRA[i]
    @inbounds dRAin[i] = β*RAout[i] - kp*RAin[i] - kdeg*(RAR[i]/(gamma+RAR[i]))*RAin[i] - mon*RAin[i].*BP[i] + moff*RABP[i] - rdeg1*RAin[i]
    @inbounds dR[i]    = Vr - rdeg2*R[i]- jalpha*RABP[i].*R[i] + jbeta*BP[i].*RAR[i]
    @inbounds dRAR[i]    = jalpha*RABP[i].*R[i] - jbeta*BP[i].*RAR[i]
    @inbounds dBP[i]     = Vbp - bpdeg1*BP[i] - mon*RAin[i].*BP[i] + moff*RABP[i] + jalpha*RABP[i].*R[i] - jbeta*BP[i].*RAR[i] + ((d.*RAR[i])./(e+RAR[i]))
    @inbounds dRABP[i]   = mon*RAin[i]*BP[i] - moff*RABP[i] - jalpha*RABP[i]*R[i] + jbeta*BP[i]*RAR[i]
  end
  nothing
end

################ Add Stochasticity

println("Solve the gradient SDE")
function ra_noise(t,u,du)
  RAout = @view u[:,:,1]
  RAin = @view u[:,:,2]
  RAR = @view u[:,:,4]
  dRAout = @view du[:,:,1]
  dRAin = @view du[:,:,2]
  dRAR = @view du[:,:,4]

  dRAout .= ϵout .* RAout
  dRAin  .= ϵin  .* RAin
  dRAR  .= ϵR  .* RAR
  nothing
end

u0 = zeros(M,N,6)
tspan2 = (0.0,500.0)
prob = SDEProblem(ra_gradient,ra_noise,u0,tspan2)

println("Solve with SOSRI")
@time noisy_ra_sol = solve(prob,SOSRI(),
                           save_everystep=false,progress_steps=10_000,
                           progress=true,abstol=1e-1,reltol=1e-2)

@time noisy_ra_sol = solve(prob,SRIW1(),
                           save_everystep=false,progress_steps=10_000,
                           progress=true,abstol=1e-5,reltol=1e-3)

@time noisy_ra_sol = solve(prob,SOSRI2(),
                           save_everystep=false,progress_steps=10_000,
                           progress=true,abstol=1e-3,reltol=1e-3)

println("Solve with EM")
# dt = 1/10000 fails
for i in 1:10
    @time noisy_ra_sol = solve(prob,EM(),dt=1/20000,
                              save_everystep=false,progress_steps=10_000,
                              progress=true)
end

# Instead try to estimate the total time

tspan3 = (0.0,1/100)
prob2 = SDEProblem(ra_gradient,ra_noise,u0,tspan3)

# LOL No.

@time noisy_ra_sol = solve(prob2,ImplicitEM(),dt=1/80000,
                          save_everystep=false,progress_steps=1,
                          progress=true)

@time noisy_ra_sol = solve(prob2,ImplicitRKMil(),dt=1/80000,
                        save_everystep=false,progress_steps=1,
                        progress=true)

################ Solve ODE

println("Solve the ODE")
tspan = (0.0,500.0)
prob_ode = ODEProblem(ra_gradient,u0,tspan)
@time ra_sol = solve(prob_ode,CVODE_BDF(linear_solver = :GMRES),save_everystep=false,
                     abstol=1e-7,reltol=1e-4)

du = zeros(N,M,6)
ra_gradient(0.0,ra_sol[end],du)
maximum(du)
