using SRKGenerator
#Standard Driver Script

##Parallelize correctly
host = gethostname()
if host == "ChrisRackBed"
  #cd("D:\\OneDrive\\Current\\projectCode\\ImprovedSRK\\Optimization")
  #numCores = 4
  #addprocs(numCores-length(workers()))
  ptx_str = "integrationWin.ptx"
  cudaCores=[1664,640]
end

if host == "crackauc.math.uci.edu"
  #numCores = 2
  #addprocs(numCores-length(workers()))
  ptx_str = "integrationLab.ptx"
  cudaCores = [2816]
end

if host == "ChrisRackTV"
  #numCores = 4
  #addprocs(numCores-length(workers()))
  ptx_str = joinpath(Pkg.dir("SRKGenerator"),"deps","integrationWin.ptx")
  cudaCores=[2560]
end

resString = srk_optimize(:GN_ISRES,1/100,Int(1e8),10000,-12,1,-1,88,5,6,gpuEnabled = true,cudaCores=cudaCores,ptx_str=ptx_str)
