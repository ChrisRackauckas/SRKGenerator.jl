using SpecialMatrices
using EllipsisNotation
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
const ϵout = 0.005
const ϵin = 0.005
const ϵR = 0.005

############ Setup Diffusion Matrices

Ax = -full(Strang(N))
Ay = -full(Strang(M))
Ax[2,1] = 2
Ax[end-1,end] = 2
Ax[end,end] = -2*(1+dx*kA)
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

  diffRA = DRA.*(RAout*Ax + Ay*RAout)
  @. dRAout  = VRA - β*RAout + kp*RAin + diffRA
  @. dRAin = β*RAout - kp*RAin - kdeg*(RAR/(gamma+RAR))*RAin - mon*RAin.*BP + moff*RABP - rdeg1*RAin
  @. dR      = Vr - rdeg2*R- jalpha*RABP.*R + jbeta*BP.*RAR
  @. dRAR    = jalpha*RABP.*R - jbeta*BP.*RAR
  @. dBP     = Vbp - bpdeg1*BP - mon*RAin.*BP + moff*RABP + jalpha*RABP.*R - jbeta*BP.*RAR + ((d.*RAR)./(e+RAR))
  @. dRABP   = mon*RAin*BP - moff*RABP - jalpha*RABP*R + jbeta*BP*RAR
  nothing
end

################ Solve ODE

println("Solve the ODE")
tspan = (0.0,5000.0)
prob = ODEProblem(ra_gradient,u0,tspan)
@time ra_sol = solve(prob,CVODE_BDF(linear_solver = :GMRES),save_everystep=false,
                     abstol=1e-7,reltol=1e-4)

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

u0 = ra_sol[end]
tspan2 = (0.0,50.0)
prob = SDEProblem(ra_gradient,ra_noise,u0,tspan2)

@time noisy_ra_sol = solve(prob,SRIW1(),
                           save_everystep=false,progress_steps=10_000,
                           progress=true,abstol=1e-5,reltol=1e-3)

@time noisy_ra_sol = solve(prob,SOSRI(),
                           save_everystep=false,progress_steps=10_000,
                           progress=true,abstol=1e-1,reltol=1e-1)

@time noisy_ra_sol = solve(prob,SOSRI2(),
                           save_everystep=false,progress_steps=10_000,
                           progress=true,abstol=1e-1,reltol=1e-2)

@time noisy_ra_sol = solve(prob,EM(),dt=1/200000,
                           save_everystep=false,progress_steps=10_000,
                           progress=true)
