using SRKGenerator

dx = 1/50
imin,imax,jmin,jmax = -12,1,-10,10
x = [-0.00575177,3.48225,-2.59676,-0.364729,1.14808,-0.783625,0.13798,-0.397231,0.668615,-0.299582,1.19174,0.107113,-0.0913619,1.58563,0.061024,-4.21259,0.605065,3.15894,-0.372408,2.31455,-1.22656,1.19005,-0.694235,-0.681141,0.97351,-0.454075,0.560372,-0.0796419,-1.82528,1.5744,0.64214,0.608841,1.6336,-2.09727,0.239234,0.224391,2.82355,-1.57304,-0.64245,-0.608839,-4.40931,4.06816,1.23893,-0.897945]
# 13.18955568047541

@time f1 = SRKGenerator.cpu_f(x,dx,imin,imax,jmin,jmax)
@time f2 = SRKGenerator.parallel_f(x,dx,imin,imax,jmin,jmax)
@time f2 = SRKGenerator.threaded_f(x,dx,imin,imax,jmin,jmax)

f1 == f2


using CUDArt

NLoptRandSeed = 0
parameter_minmax=5
max_eval=Int(1e8)
initCon = ones(44)
tol = 1e-2
ftol = 1e-15
tol2 = 1e-5
counterSteps=Int(1e5)
counterSteps2=Int(1e6)
initStepSize=[]
gpuEnabled=true
ptx_str  = joinpath(Pkg.dir("SRKGenerator"),"deps","integrationWin.ptx")
cudaCores = 1664
initStepSize2=1e-6
outfile=""
constrain_c = true

N = 26
N2= 16
M = 44
count = [0]

## Script Start
x_L = -parameter_minmax*ones(M)
x_U =  parameter_minmax*ones(M)
g_L = -tol*ones(N)
g_U =  tol*ones(N)
sizei = length(imin:dx:imax)
sizej = length(jmin:dx:jmax)
totArea = (imax-imin)*(jmax-jmin)/(sizei*sizej)
if gpuEnabled
  iarr = convert(Vector{Float32},collect(imin:dx:imax))
  jarr = convert(Vector{Float32},collect(jmin:dx:jmax))
  numCards = length(cudaCores)
  totCores = sum(cudaCores)
  equalDiv = sizei*sizej÷totCores + 1
  portions = equalDiv*cudaCores
  if numCards > 1
    startIdx = [0 cumsum(portions)[1:end-1]]
  else
    startIdx = [0]
  end
  g_iarr = Vector{CUDArt.CudaArray{Float32,1}}(numCards)
  g_jarr = Vector{CUDArt.CudaArray{Float32,1}}(numCards)
  g_tmp = Vector{CUDArt.CudaArray{Int32,1}}(numCards)
  g_coefs = Vector{CUDArt.CudaArray{Float32,1}}(numCards)
  ans = Vector{Int32}(numCards)
  integrationFuncs = Vector{CUDArt.CuFunction}(numCards)
  CUDArt.init(devices(dev->true))
  devlist = devices(dev->true)
  streams = [(device(dev); Stream()) for dev in devlist]
  for i = 1:numCards
    device(i-1)
    g_iarr[i]  = CudaArray(iarr)
    g_jarr[i]  = CudaArray(jarr)
    g_tmp[i]   = CudaArray(Int32,cudaCores[i])
    md = CuModule(ptx_str,false)
    integrationFuncs[i] = CuFunction(md,"integration")
  end
end

eV = [1;1;1;1]

eval_f = (x,grad) -> SRKGenerator.f_maker(x,ans,integrationFuncs,cudaCores,numCards,g_coefs,g_iarr,g_jarr,sizei,sizej,equalDiv,startIdx,g_tmp,totArea,counterSteps,counterSteps2,outfile,gpuEnabled,count)
eval_g = (tmp,x,grad) -> SRKGenerator.g_maker(x,tmp,eV,counterSteps,counterSteps2,outfile,count)
eval_g_ineq = (tmp,x,grad) -> SRKGenerator.g_ineq_maker(x,tmp,eV,counterSteps,counterSteps2,outfile,count)

f3 = eval_f(x,0)

### Now test some gradients

f = (x) -> SRKGenerator.cpu_f(x,dx,imin,imax,jmin,jmax)

using ForwardDiff

ForwardDiff.gradient(f,x)

using Calculus

Calculus.gradient(f,x)

gpu_f = (x) -> eval_f(x,0)

Calculus.gradient(gpu_f,x)
