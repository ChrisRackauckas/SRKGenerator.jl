function cpu_f(x,dx,imin,imax,jmin,jmax)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  coefs,powz,poww = getCoefficients(A0,A1,B0,B1,α,β1,β2,β3,β4)
  sizei = length(imin:dx:imax)
  sizej = length(jmin:dx:jmax)
  totArea = (imax-imin)*(jmax-jmin)/(sizei*sizej)
  tmp = 0
  for i = imin:dx:imax
    isq2 = i*i; isq3 = i*isq2; isq4 = isq2*isq2; isq5 = i*isq4
    isq6 = isq4*isq2; isq7 = i*isq6; isq8 = isq4*isq4
    for j=jmin:dx:jmax
      jsq2 = j*j; jsq3= j*jsq2; jsq4 = jsq2*jsq2;
      jsq5 = j*jsq4; jsq6 = jsq2*jsq4; jsq7 = j*jsq6; jsq8 = jsq4*jsq4
      @inbounds tmp += abs(coefs[1]*(jsq2) + coefs[2]*(jsq3) + coefs[3]*(jsq4) + coefs[4]*(jsq5) + coefs[5]*jsq6 + coefs[6]*jsq7 + coefs[7]*jsq8 + coefs[8]*(i) + coefs[9]*(i)*(jsq2) +
      coefs[10]*i*jsq3 + coefs[11]*(i)*(jsq4) + coefs[12]*i*jsq5 + coefs[13]*(i)*(jsq6) + coefs[14]*i*jsq7 + coefs[15]*(isq2) + coefs[16]*(isq2)*(jsq2) + coefs[17]*isq2*jsq3 +
      coefs[18]*(isq2)*(jsq4) + coefs[19]*isq2*jsq5 + coefs[20]*(isq2)*(jsq6) + coefs[21]*(isq3) + coefs[22]*(isq3)*(jsq2) + coefs[23]*isq3*jsq3 + coefs[24]*(isq3)*(jsq4) + coefs[25]*isq3*jsq5 +
      coefs[26]*(isq4) + coefs[27]*(isq4)*(jsq2) + coefs[28]*isq4*jsq3 + coefs[29]*(isq4)*(jsq4) + coefs[30]*(isq5) + coefs[31]*(isq5)*(jsq2) + coefs[32]*isq5*jsq3+ coefs[33]*(isq6) +
      coefs[34]*(isq6)*(jsq2) + coefs[35]*(isq7) + coefs[36]*(isq8))<1
    end
  end
  res = tmp*totArea
end

function threaded_f(x,dx,imin,imax,jmin,jmax)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  coefs,powz,poww = getCoefficients(A0,A1,B0,B1,α,β1,β2,β3,β4)
  sizei = length(imin:dx:imax)
  sizej = length(jmin:dx:jmax)
  totArea = (imax-imin)*(jmax-jmin)/(sizei*sizej)
  tmps = Vector{Int}(Threads.nthreads())
  tmps .= 0
  Threads.@threads for i = imin:dx:imax
    isq2 = i*i; isq3 = i*isq2; isq4 = isq2*isq2; isq5 = i*isq4
    isq6 = isq4*isq2; isq7 = i*isq6; isq8 = isq4*isq4
    for j=jmin:dx:jmax
      jsq2 = j*j; jsq3= j*jsq2; jsq4 = jsq2*jsq2;
      jsq5 = j*jsq4; jsq6 = jsq2*jsq4; jsq7 = j*jsq6; jsq8 = jsq4*jsq4
      @inbounds tmps[Threads.threadid()] += abs(coefs[1]*(jsq2) + coefs[2]*(jsq3) + coefs[3]*(jsq4) + coefs[4]*(jsq5) + coefs[5]*jsq6 + coefs[6]*jsq7 + coefs[7]*jsq8 + coefs[8]*(i) + coefs[9]*(i)*(jsq2) +
      coefs[10]*i*jsq3 + coefs[11]*(i)*(jsq4) + coefs[12]*i*jsq5 + coefs[13]*(i)*(jsq6) + coefs[14]*i*jsq7 + coefs[15]*(isq2) + coefs[16]*(isq2)*(jsq2) + coefs[17]*isq2*jsq3 +
      coefs[18]*(isq2)*(jsq4) + coefs[19]*isq2*jsq5 + coefs[20]*(isq2)*(jsq6) + coefs[21]*(isq3) + coefs[22]*(isq3)*(jsq2) + coefs[23]*isq3*jsq3 + coefs[24]*(isq3)*(jsq4) + coefs[25]*isq3*jsq5 +
      coefs[26]*(isq4) + coefs[27]*(isq4)*(jsq2) + coefs[28]*isq4*jsq3 + coefs[29]*(isq4)*(jsq4) + coefs[30]*(isq5) + coefs[31]*(isq5)*(jsq2) + coefs[32]*isq5*jsq3+ coefs[33]*(isq6) +
      coefs[34]*(isq6)*(jsq2) + coefs[35]*(isq7) + coefs[36]*(isq8))<1
    end
  end
  tmp = sum(tmps)
  res = tmp*totArea
end

function parallel_f(x,dx,imin,imax,jmin,jmax)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  coefs,powz,poww = getCoefficients(A0,A1,B0,B1,α,β1,β2,β3,β4)
  sizei = length(imin:dx:imax)
  sizej = length(jmin:dx:jmax)
  totArea = (imax-imin)*(jmax-jmin)/(sizei*sizej)
  ans = @sync @parallel (+) for i = imin:dx:imax
    tmp = 0
    isq2 = i*i; isq3 = i*isq2; isq4 = isq2*isq2; isq5 = i*isq4
    isq6 = isq4*isq2; isq7 = i*isq6; isq8 = isq4*isq4
    for j=jmin:dx:jmax
      jsq2 = j*j; jsq3= j*jsq2; jsq4 = jsq2*jsq2;
      jsq5 = j*jsq4; jsq6 = jsq2*jsq4; jsq7 = j*jsq6; jsq8 = jsq4*jsq4
      @inbounds tmp += abs(coefs[1]*(jsq2) + coefs[2]*(jsq3) + coefs[3]*(jsq4) + coefs[4]*(jsq5) + coefs[5]*jsq6 + coefs[6]*jsq7 + coefs[7]*jsq8 + coefs[8]*(i) + coefs[9]*(i)*(jsq2) +
      coefs[10]*i*jsq3 + coefs[11]*(i)*(jsq4) + coefs[12]*i*jsq5 + coefs[13]*(i)*(jsq6) + coefs[14]*i*jsq7 + coefs[15]*(isq2) + coefs[16]*(isq2)*(jsq2) + coefs[17]*isq2*jsq3 +
      coefs[18]*(isq2)*(jsq4) + coefs[19]*isq2*jsq5 + coefs[20]*(isq2)*(jsq6) + coefs[21]*(isq3) + coefs[22]*(isq3)*(jsq2) + coefs[23]*isq3*jsq3 + coefs[24]*(isq3)*(jsq4) + coefs[25]*isq3*jsq5 +
      coefs[26]*(isq4) + coefs[27]*(isq4)*(jsq2) + coefs[28]*isq4*jsq3 + coefs[29]*(isq4)*(jsq4) + coefs[30]*(isq5) + coefs[31]*(isq5)*(jsq2) + coefs[32]*isq5*jsq3+ coefs[33]*(isq6) +
      coefs[34]*(isq6)*(jsq2) + coefs[35]*(isq7) + coefs[36]*(isq8))<1
    end
    tmp
  end
  res = sum(ans)*totArea
end

function f_maker(x,ans,coefs,powz,poww,integrationFuncs,cudaCores,numCards,g_coefs,g_iarr,g_jarr,sizei,sizej,equalDiv,startIdx,g_tmp,totArea,counterSteps,counterSteps2,outfile,gpuEnabled,count)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  getCoefficients!(coefs,powz,poww,A0,A1,B0,B1,α,β1,β2,β3,β4)
  if gpuEnabled # Commented out parts for multiple graphics cards
  #  @sync begin
      for i = 1:numCards
      #  @async begin
          dev = device(i-1)
          g_coefs[i] = CudaArray(coefs)
          launch(integrationFuncs[i], cudaCores[i], 64, (g_coefs[i], g_iarr[i], g_jarr[i], sizei,sizej,equalDiv,startIdx[i],g_tmp[i],))#,stream=streams[i]);
          #wait(streams[i])
          ans[i] = sum(to_host(g_tmp[i]))
      #  end
      end
    # end
    res = sum(ans)*totArea
  else
    ans = @sync @parallel (+) for i = imin:dx:imax
      tmp = 0
      isq2 = i*i; isq3 = i*isq2; isq4 = isq2*isq2; isq5 = i*isq4
      isq6 = isq4*isq2; isq7 = i*isq6; isq8 = isq4*isq4
      for j=jmin:dx:jmax
        jsq2 = j*j; jsq3= j*jsq2; jsq4 = jsq2*jsq2;
        jsq5 = j*jsq4; jsq6 = jsq2*jsq4; jsq7 = j*jsq6; jsq8 = jsq4*jsq4
        @inbounds tmp += abs(coefs[1]*(jsq2) + coefs[2]*(jsq3) + coefs[3]*(jsq4) + coefs[4]*(jsq5) + coefs[5]*jsq6 + coefs[6]*jsq7 + coefs[7]*jsq8 + coefs[8]*(i) + coefs[9]*(i)*(jsq2) +
        coefs[10]*i*jsq3 + coefs[11]*(i)*(jsq4) + coefs[12]*i*jsq5 + coefs[13]*(i)*(jsq6) + coefs[14]*i*jsq7 + coefs[15]*(isq2) + coefs[16]*(isq2)*(jsq2) + coefs[17]*isq2*jsq3 +
        coefs[18]*(isq2)*(jsq4) + coefs[19]*isq2*jsq5 + coefs[20]*(isq2)*(jsq6) + coefs[21]*(isq3) + coefs[22]*(isq3)*(jsq2) + coefs[23]*isq3*jsq3 + coefs[24]*(isq3)*(jsq4) + coefs[25]*isq3*jsq5 +
        coefs[26]*(isq4) + coefs[27]*(isq4)*(jsq2) + coefs[28]*isq4*jsq3 + coefs[29]*(isq4)*(jsq4) + coefs[30]*(isq5) + coefs[31]*(isq5)*(jsq2) + coefs[32]*isq5*jsq3+ coefs[33]*(isq6) +
        coefs[34]*(isq6)*(jsq2) + coefs[35]*(isq7) + coefs[36]*(isq8))<1
      end
      tmp
    end
  end
  count[1]=count[1]+1
  if mod(count[1],counterSteps)==0
    prints = count÷counterSteps
    println("f_$(prints[1])($x)")
    println("$res")
    flush(STDOUT)
  end
  if mod(count[1],counterSteps2)==0
    mathx = translateToMathematica(x)
    println(mathx)
    flush(STDOUT)
    if outfile != ""
      prints = count÷counterSteps
      julString = printForJulia(x)
      saveStr = """

      ------------------------------------------------------------------------
      f_$(prints[1])($x)

      $res

      $mathx

      $julString

      """
      write(outfile,saveStr)
      flush(outfile)
    end
  end
  return res
end

function g_maker(x,tmp,eV,counterSteps,counterSteps2,outfile,count)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  tmp[1] = sum(α)-1
  tmp[2] = sum(β1)-1
  tmp[3] = sum(β2)
  tmp[4] = sum(β3)
  tmp[5] = sum(β4)
  tmp[6] = (β1'*B1*eV)[1]
  tmp[7] = (transpose(β2)*B1*eV-1)[1]
  tmp[9] = (transpose(β3)*B1*eV)[1]
  tmp[10]= (transpose(β4)*B1*eV)[1]
  tmp[11]= (transpose(α)*A0*eV-1/2)[1]
  tmp[12]= (transpose(α)*B0*eV-1)[1]
  tmp[13]= (transpose(α)*(B0*eV).^2-3/2)[1]
  tmp[14]= (transpose(β1)*A1*eV-1)[1]
  tmp[15]= (transpose(β2)*A1*eV)[1]
  tmp[16]= (transpose(β3)*A1*eV+1)[1]
  tmp[17]= (transpose(β4)*A1*eV)[1]
  tmp[18]= (transpose(β1)*(B1*eV).^2-1)[1]
  tmp[19]= (transpose(β2)*(B1*eV).^2)[1]
  tmp[20]= (transpose(β3)*(B1*eV).^2+1)[1]
  tmp[21]= (transpose(β4)*(B1*eV).^2-2)[1]
  tmp[22]=(transpose(β1)*(B1*(B1*eV))-0)[1]
  tmp[23]=(transpose(β2)*(B1*(B1*eV))-0)[1]
  tmp[24]=(transpose(β3)*(B1*(B1*eV))-0)[1]
  tmp[25]=(transpose(β4)*(B1*(B1*eV))-1)[1]
  tmp[26]=((1//2)*transpose(β1)*(A1*(B0*eV))+(1//3)*transpose(β3)*(A1*(B0*eV)))[1]
  maxErr = maximum(tmp)
  if maxErr > 1e10
    error("Diverged")
  end
  if mod(count[1],counterSteps)==0
    println("Max g: $maxErr")
    flush(STDOUT)
  end
  if mod(count[1],counterSteps2)==0
    if outfile != ""
      maxErr = maximum(tmp)
      write(outfile,"Max g: $maxErr")
      flush(outfile)
    end
  end
  nothing
end

function g_ineq_maker(x,tmp,eV,counterSteps,counterSteps2,outfile,count)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  c0 = A0*eV
  c1 = A1*eV
  tmp[1:4] .= c0 .- 1
  tmp[5:8] .= -(c0)
  tmp[9:12] .= c1 .- 1
  tmp[13:16] .= -(c1)
  nothing
end

#Dg = ForwardDiff.jacobian(g)
#∇f = ForwardDiff.gradient(f)

function eval_g(result,x,grad)
  #=
  if length(grad)>0
    grad[:]=Dg(x)
  end
  =#
  g(x,result)
end


#function eval_f(x,grad)
  #=
  if length(grad)>0
    grad[:]=∇f(x)
  end
  =#
#  return(f(x))
#end
