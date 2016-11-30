function f(x)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  coefs,powz,poww = getCoefficients(A0,A1,B0,B1,α,β1,β2,β3,β4)
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
      @simd for j=jmin:dx:jmax
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
  res = sum(ans)*totArea
  count+=1
  if mod(count,counterSteps)==0
    prints = count÷counterSteps
    println("f_$prints($x)")
    println("$res")
    #println(whos())
    #println(whos(Base))
    #=
    for name in setdiff(names(Base, true), names(Base))
         try
             s = Base.summarysize(getfield(Base, name))
             if s > 10000000
                 println(name, " ", s)
             end
         end
     end
    =#
    flush(STDOUT)
  end
  if mod(count,counterSteps2)==0
    mathx = translateToMathematica(x)
    println(mathx)
    flush(STDOUT)
    if outfile != ""
      prints = count÷counterSteps
      julString = printForJulia(x)
      saveStr = """

      ------------------------------------------------------------------------
      f_$prints($x)

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

function g(x)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  eV = [1;1;1;1]
  ans = zeros(eltype(x),26)
  ans[1] = sum(α)-1
  ans[2] = sum(β1)-1
  ans[3] = sum(β2)
  ans[4] = sum(β3)
  ans[5] = sum(β4)
  ans[6] = (β1'*B1*eV)[1]
  ans[7] = (transpose(β2)*B1*eV-1)[1]
  ans[9] = (transpose(β3)*B1*eV)[1]
  ans[10]= (transpose(β4)*B1*eV)[1]
  ans[11]= (transpose(α)*A0*eV-1/2)[1]
  ans[12]= (transpose(α)*B0*eV-1)[1]
  ans[13]= (transpose(α)*(B0*eV).^2-3/2)[1]
  ans[14]= (transpose(β1)*A1*eV-1)[1]
  ans[15]= (transpose(β2)*A1*eV)[1]
  ans[16]= (transpose(β3)*A1*eV+1)[1]
  ans[17]= (transpose(β4)*A1*eV)[1]
  ans[18]= (transpose(β1)*(B1*eV).^2-1)[1]
  ans[19]= (transpose(β2)*(B1*eV).^2)[1]
  ans[20]= (transpose(β3)*(B1*eV).^2+1)[1]
  ans[21]= (transpose(β4)*(B1*eV).^2-2)[1]
  ans[22]=(transpose(β1)*(B1*(B1*eV))-0)[1]
  ans[23]=(transpose(β2)*(B1*(B1*eV))-0)[1]
  ans[24]=(transpose(β3)*(B1*(B1*eV))-0)[1]
  ans[25]=(transpose(β4)*(B1*(B1*eV))-1)[1]
  ans[26]=((1//2)*transpose(β1)*(A1*(B0*eV))+(1//3)*transpose(β3)*(A1*(B0*eV)))[1]
  if mod(count,counterSteps)==0
    maxErr = maximum(ans)
    println("Max g: $maxErr")
    flush(STDOUT)
  end
  if mod(count,counterSteps2)==0
    if outfile != ""
      maxErr = maximum(ans)
      write(outfile,"Max g: $maxErr")
      flush(outfile)
    end
  end
  return ans
end

#Dg = ForwardDiff.jacobian(g)
#∇f = ForwardDiff.gradient(f)

function eval_g(result,x,grad)
  #=
  if length(grad)>0
    grad[:]=Dg(x)
  end
  =#
  result[:]=g(x)
end

function eval_f(x,grad)
  #=
  if length(grad)>0
    grad[:]=∇f(x)
  end
  =#
  return(f(x))
end