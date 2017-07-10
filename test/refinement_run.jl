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

initCon =  [0.13803991182809405,0.5764629771061757,0.4235370228938241,0.45988142444172125,0.8084859191118128,-0.2683673489095204,0.45482884244918836,0.7519704601202507,0.24802953987974916,0.6961685917464052,0.3460790565637937,-0.04224764831019894,0.08513289389320751,1.0272930820405226,0.4565894293733811,1.7266402546369712,-0.46040416640520754,-0.9623742513556953,0.6744100174205603,-0.07377591084107533,-0.4985455217486592,-0.5623853411970345,0.023606441329882558,-0.8953332413687689,-0.14925077940925455,0.7532260416025327,0.6911225838401098,-0.2950978480741201,-0.4550650633914512,0.8347196075716,0.3790963153683881,0.24124914164997122,-0.5014671170389827,0.9198342527636846,-0.2556662881012024,-0.16270085567402398,1.4550650304633204,-0.8347196116875014,-0.37909625532888075,-0.24124916721822443,-0.4998591615011572,0.9168848205568053,-1.4115003346202473,0.994474673400943]
x_L = initCon .- 0.01
x_U = initCon .+ 0.01
resString = srk_optimize(:LN_AUGLAG_EQ,1/100,0,imin,1,-jmax,jmax,
                         rand_minmax = rand(1:3),tol=1e-15,
                         initCon = initCon,x_L=x_L,x_U=x_U,
                         NLoptRandSeed=rand(1:Int(1e8)),gpuEnabled = true,
                         cudaCores=cudaCores,ptx_str=ptx_str)
