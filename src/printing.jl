function translateToMathematica(x)
  rhs = """Set[{A021,A031,A032,A041,A042,A043,A121,A131,A132,A141,A142,A143,B021,B031,B032,B041,B042,B043,B121,B131,B132,B141,B142,B143,\\[Alpha]1,\\[Alpha]2,\\[Alpha]3,\\[Alpha]4,\\[Beta]11,\\[Beta]12,\\[Beta]13,\\[Beta]14,\\[Beta]21,\\[Beta]22,\\[Beta]23,\\[Beta]24,\\[Beta]31,\\[Beta]32,\\[Beta]33,\\[Beta]34,\\[Beta]41,\\[Beta]42,\\[Beta]43,\\[Beta]44}"""

  lhs = @sprintf("{%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s}]",x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],x[30],x[31],x[32],x[33],x[34],x[35],x[36],x[37],x[38],x[39],x[40],x[41],x[42],x[43],x[44])
  println(lhs)
  mathstr = "$rhs,$lhs\n"
  return mathstr
end

function MathematicaInitSet(x)
  lhs = @sprintf("{A021,%s},{A031,%s},{A032,%s},{A041,%s},{A042,%s},{A043,%s},{A121,%s},{A131,%s},{A132,%s},{A141,%s},{A142,%s},{A143,%s},{B021,%s},{B031,%s},{B032,%s},{B041,%s},{B042,%s},{B043,%s},{B121,%s},{B131,%s},{B132,%s},{B141,%s},{B142,%s},{B143,%s},{\\[Alpha]1,%s},{\\[Alpha]2,%s},{\\[Alpha]3,%s},{\\[Alpha]4,%s},{\\[Beta]11,%s},{\\[Beta]12,%s},{\\[Beta]13,%s},{\\[Beta]14,%s},{\\[Beta]21,%s},{\\[Beta]22,%s},{\\[Beta]23,%s},{\\[Beta]24,%s},{\\[Beta]31,%s},{\\[Beta]32,%s},{\\[Beta]33,%s},{\\[Beta]34,%s},{\\[Beta]41,%s},{\\[Beta]42,%s},{\\[Beta]43,%s},{\\[Beta]44,%s}\n",x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],x[30],x[31],x[32],x[33],x[34],x[35],x[36],x[37],x[38],x[39],x[40],x[41],x[42],x[43],x[44])
  println(lhs)
  return lhs
end

function MathematicaReplaceSet(x)
  lhs = @sprintf("{A021->%s,A031->%s,A032->%s,A041->%s,A042->%s,A043->%s,A121->%s,A131->%s,A132->%s,A141->%s,A142->%s,A143->%s,B021->%s,B031->%s,B032->%s,B041->%s,B042->%s,B043->%s,B121->%s,B131->%s,B132->%s,B141->%s,B142->%s,B143->%s,\\[Alpha]1->%s,\\[Alpha]2->%s,\\[Alpha]3->%s,\\[Alpha]4->%s,\\[Beta]11->%s,\\[Beta]12->%s,\\[Beta]13->%s,\\[Beta]14->%s,\\[Beta]21->%s,\\[Beta]22->%s,\\[Beta]23->%s,\\[Beta]24->%s,\\[Beta]31->%s,\\[Beta]32->%s,\\[Beta]33->%s,\\[Beta]34->%s,\\[Beta]41->%s,\\[Beta]42->%s,\\[Beta]43->%s,\\[Beta]44->%s}\n",x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],x[30],x[31],x[32],x[33],x[34],x[35],x[36],x[37],x[38],x[39],x[40],x[41],x[42],x[43],x[44])
  println(lhs)
  return lhs
end


function printForJulia(x)
  A0,A1,B0,B1,α,β1,β2,β3,β4 = translate(x)
  c0 = A0*eV
  c1 = A1*eV
  resString = """
  x = [$(join(string.(x), ","))]
  c0 = [$(join(string.(c0), ","))]
  c1 = [$(join(string.(c1), ","))]
  A0 = [$(join(string.(A0), ","))]
  A1 = [$(join(string.(A1), ","))]
  B0 = [$(join(string.(B0), ","))]
  B1 = [$(join(string.(B1), ","))]
  α  = [$(join(string.(α), ","))]
  β1 = [$(join(string.(β1), ","))]
  β2 = [$(join(string.(β2), ","))]
  β3 = [$(join(string.(β3), ","))]
  β4 = [$(join(string.(β4), ","))]
  """
  return resString
end
