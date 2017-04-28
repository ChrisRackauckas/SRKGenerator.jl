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
  cudaCores = [parse(Int,ARGS[1])]
  imin = parse(Int,ARGS[2])
  jmax = parse(Int,ARGS[3]); jmin = -jmax


end

resString = srk_optimize(:LN_AUGLAG_EQ,1/100,0,imin,1,-jmax,jmax,
                         rand_minmax = rand(1:3),
                         NLoptRandSeed=rand(1:Int(1e8)),gpuEnabled = true,
                         cudaCores=cudaCores,ptx_str=ptx_str)
