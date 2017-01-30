function srk_optimize(alg,dx,mev,populationSize,imin,imax,jmin,jmax,len,NLoptRandSeed;
                    initCon = ones(44),tol = 1e-2,ftol = 1e-15,tol2 = 1e-5,
                    counterSteps=Int(1e5),counterSteps2=Int(1e6),
                    initStepSize=[],gpuEnabled=true,ptx_str  = "integration.ptx",
                    cudaCores = 1664,initStepSize2=1e-6,outfile="")
  ##Parameters
  N = 26
  M = 44
  count = [0]
  timeNow = now()
  timeNow = Dates.format(timeNow, "y-m-d-HH-MM-SS")
  outfile = open(joinpath(Pkg.dir("SRKGenerator"),"output","Opti$timeNow.txt"), "w")

  ## Script Start
  x_L = -len*ones(M)
  x_U =  len*ones(M)
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

  eval_f = (x,grad) -> f_maker(x,ans,integrationFuncs,cudaCores,numCards,g_coefs,g_iarr,g_jarr,sizei,sizej,equalDiv,startIdx,g_tmp,totArea,counterSteps,counterSteps2,outfile,gpuEnabled,count)
  eval_g = (tmp,x,grad) -> g_maker(x,tmp,eV,counterSteps,counterSteps2,outfile,count)

  opt = Opt(alg,M) #:LD_SLSQP, :LN_COBYLA (semi), :GN_ISRES support equality constraints
  lower_bounds!(opt,x_L)
  upper_bounds!(opt,x_U)
  max_objective!(opt, eval_f)
  maxeval!(opt::Opt, mev::Integer)
  equality_constraint!(opt,eval_g,tol*ones(N))
  ftol_abs!(opt,ftol)
  xtol_abs!(opt,ftol)
  if initStepSize != []
    initial_step!(opt,initStepSize)
  end
  population!(opt,populationSize)
  NLopt.srand(NLoptRandSeed)

  minf,minx,ret = optimize(opt,initCon)
  resString = """


  -----------------Final Result------------------
  Options: alg=$alg,dx=$dx,mev=$mev,popSize=$populationSize
  imin=$imin,jmin=$jmin,imax=$imax,jmax=$jmax,len=$len,randSeed=$NLoptRandSeed
  """
  println(resString)
  if gpuEnabled
    CUDArt.close(devices(dev->true))
  end
  write(outfile,resString)
  flush(outfile)

  return resString
end
