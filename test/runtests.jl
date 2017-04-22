using SRKGenerator
#Standard Driver Script

if isempty(ARGS)

  ## Called from test. Use a default setup
  imin = -12
  jmax = 1

  ## Parallelize correctly
  host = gethostname()
  if host == "ChrisRackBed"
    #numCores = 4
    #addprocs(numCores-length(workers()))
    ptx_str = joinpath(Pkg.dir("SRKGenerator"),"deps","integrationWin.ptx")
    cudaCores=[1664,640]
  end

  if host == "crackauc2"
    #numCores = 16
    #addprocs(numCores-length(workers()))
    ptx_str = joinpath(Pkg.dir("SRKGenerator"),"deps","integrationLinux.ptx")
    cudaCores = [2816]
  end

  if host == "ChrisRackTV"
    #numCores = 4
    #addprocs(numCores-length(workers()))
    ptx_str = joinpath(Pkg.dir("SRKGenerator"),"deps","integrationWin.ptx")
    cudaCores=[2560]
  end

else

  ## Called from commandline. Assume Linux.
  ptx_str = joinpath(Pkg.dir("SRKGenerator"),"deps","integrationLinux.ptx")
  ## Arg Form:
  cudaCores = ARGS[1]
  imin = ARGS[2]
  jmax = ARGS[3]; jmin = -jmax


end

resString = srk_optimize(:LN_AUGLAG_EQ,1/100,100000,-imin,1,-jmax,jmax,
                         NLoptRandSeed=rand(1:Int(1e8)),gpuEnabled = true,
                         cudaCores=cudaCores,ptx_str=ptx_str)
