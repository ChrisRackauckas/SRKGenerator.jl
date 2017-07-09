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

initCon =  [0.07519788885593696,-0.06949643927329126,1.0694964392681514,-0.1522251222812998,1.0998266313011555,0.05239849098014294,0.4432893980421528,0.9925661306103539,0.007433869389646297,0.7875149068373734,0.25016319442776613,-0.03767810126513897,0.1088062121063868,0.4577017766471385,0.3547124897943134,0.6269734317481137,1.262208306979323,-0.39725204936894987,-0.6657998109848107,-1.08367353263802,0.4343087675234363,3.712478464631237,-0.158810046896599,-1.7033531774688744,-0.7413336106671035,1.3422694398238284,-0.38053155726800675,0.7795957286734092,0.15796207948071708,-0.28374182014598404,0.9089014017198427,0.21687833951738097,0.9413242737406766,-1.6908682348962363,0.6051464035225309,0.1443975579842588,0.8420379278988441,0.2837418089595506,-0.9089013961726579,-0.21687833942011436,-1.5676321255689545,2.8158833702962753,-1.6740034356296143,0.42575219115133534]
resString = srk_optimize(:LN_COBYLA,1/100,0,imin,1,-jmax,jmax,
                         rand_minmax = rand(1:3),
                         initCon = initCon,
                         NLoptRandSeed=rand(1:Int(1e8)),gpuEnabled = true,
                         cudaCores=cudaCores,ptx_str=ptx_str)
