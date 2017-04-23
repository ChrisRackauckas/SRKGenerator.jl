function srk_optimize(alg,dx,pop_size,imin,imax,jmin,jmax;
                    NLoptRandSeed = 0,parameter_minmax=5,max_eval=Int(1e8),
                    initCon = ones(44),tol = 1e-2,ftol = 1e-15,tol2 = 1e-5,
                    counterSteps=Int(1e5),counterSteps2=Int(1e6),
                    initStepSize=[],gpuEnabled=true,ptx_str  = "integration.ptx",
                    cudaCores = 1664,initStepSize2=1e-6,outfile="",constrain_c = true)
  ##Parameters
  N = 26
  N2= 16
  M = 44
  count = [0]
  timeNow = now()
  timeNow = Dates.format(timeNow, "y-m-d-HH-MM-SS")
  outfile = open(joinpath(Pkg.dir("SRKGenerator"),"output","Opti$timeNow.txt"), "w")

  ## Script Start
  coef_ans = Vector{Float32}(36)
  powz = Vector{Int8}(36)
  poww = Vector{Int8}(36)
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

  eval_f = (x,grad) -> f_maker(x,ans,coef_ans,powz,poww,integrationFuncs,cudaCores,numCards,g_coefs,g_iarr,g_jarr,sizei,sizej,equalDiv,startIdx,g_tmp,totArea,counterSteps,counterSteps2,outfile,gpuEnabled,count)
  eval_g = (tmp,x,grad) -> g_maker(x,tmp,eV,counterSteps,counterSteps2,outfile,count)
  eval_g_ineq = (tmp,x,grad) -> g_ineq_maker(x,tmp,eV,counterSteps,counterSteps2,outfile,count)
  opt = Opt(alg,M) #:LD_SLSQP, :LN_COBYLA (semi), :GN_ISRES support equality constraints
  lower_bounds!(opt,x_L)
  upper_bounds!(opt,x_U)
  max_objective!(opt, eval_f)
  maxeval!(opt::Opt, max_eval::Integer)
  equality_constraint!(opt,eval_g,tol*ones(N))
  if constrain_c
    inequality_constraint!(opt,eval_g_ineq,tol*ones(N2))
  end
  ftol_abs!(opt,ftol)
  xtol_abs!(opt,ftol)
  if initStepSize != []
    initial_step!(opt,initStepSize)
  end
  population!(opt,pop_size)
  NLopt.srand(NLoptRandSeed)

  setupString = """

  -----------------Setup------------------
  Options: alg=$alg,dx=$dx,max_eval=$max_eval,popSize=$pop_size
  imin=$imin,jmin=$jmin,imax=$imax,jmax=$jmax,
  parameter_minmax=$parameter_minmax,randSeed=$NLoptRandSeed
  """
  println(setupString)
  write(outfile,setupString)
  flush(outfile)

  maxf,maxx,ret = optimize(opt,initCon)

  tmp2 = Vector{Float32}(N)
  eval_g(tmp2,maxx,0)
  maxErr = maximum(tmp2)
  mathx = translateToMathematica(maxx)
  println(mathx)
  flush(STDOUT)
  julString = printForJulia(maxx)
  resString = """

  -----------------Final Result------------------
  Options: alg=$alg,dx=$dx,max_eval=$max_eval,popSize=$pop_size
  imin=$imin,jmin=$jmin,imax=$imax,jmax=$jmax,parameter_minmax=$parameter_minmax,randSeed=$NLoptRandSeed

  $(count[1]) steps

  Completed with maxf : $maxf

  Constraints at $maxErr

  $julString

  $mathx

  """
  println(resString)
  if gpuEnabled
    CUDArt.close(devices(dev->true))
  end
  write(outfile,resString)
  flush(outfile)

  return resString
end
