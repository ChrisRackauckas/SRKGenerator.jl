using SRKGenerator
#Standard Driver Script

##Parallelize correctly
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

resString = srk_optimize(:LN_AUGLAG_EQ,1/100,100000,-12,1,-1,1,
                         NLoptRandSeed=6,gpuEnabled = true,
                         cudaCores=cudaCores,ptx_str=ptx_str)
