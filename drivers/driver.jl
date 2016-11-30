#Standard Driver Script

##Parallelize correctly
host = gethostname()
if host == "ChrisRackBed"
  cd("D:\\OneDrive\\Current\\projectCode\\ImprovedSRK\\Optimization")
  numCores = 4
  addprocs(numCores-length(workers()))
  ptxStr = "integrationWin.ptx"
  cudaCores=[1664,640]
end

if host == "crackauc.math.uci.edu"
  numCores = 2
  addprocs(numCores-length(workers()))
  ptxStr = "integrationLab.ptx"
  cudaCores = [2816]
end

if host == "ChrisRackTV"
  numCores = 4
  addprocs(numCores-length(workers()))
  ptxStr = "integrationWin.ptx"
  cudaCores=[1664]
end

## Else HPC, spawn via julia -p n
## numcores - 1 processes launched.
using SRKGenerator
timeNow = now()
timeNow = Dates.format(timeNow, "y-m-d-HH-MM-SS")
resString = SRKoptimize(:GN_ISRES,1/200,Int(5e7),4000,-2,1,-4,4,5,30,mev2=10,gpuEnabled = true,cudaCores=cudaCores,ptxStr=ptxStr)
CUDArt.close(devices(dev->true))
outfile = open("output/Opti$timeNow.txt", "w")
write(outfile,resString)
flush(outfile)
