using NLopt, ForwardDiff

x0 = zeros(10)
gcfg = ForwardDiff.GradientConfig(nothing,x0)
g!(grad, x) = ForwardDiff.gradient!(grad, integral_function, x, gcfg)
lbs  = -5ones(10)
lbs[6] = 0.0
ubs = 20ones(10)
ubs[6] = 1.0

function integral_function(x)
  a2, a3, b13, b22, b23, A021, A031, B021, B031, B032 = x
  c01 = 0
  c03 = 1
  c11 = 0
  c13 = 1
  A032 = 1 - A031
  b12 = (b22*(1 - c11 + b13*c11 - b13*c13))/(-1 + b23*c11 - b23*c13)
  b11 = 1 - b12 - b13
  b21 = - b22 - b23
  c12 = (-1 - b21*c11 - b23*c13)/b22
  a1 = 1 - a2 - a3
  c02 = A021
  G(A021,A031,A032,a1,a2,a3)
end

iter = 0
g_res = zeros(3)

function f(x::Vector, grad::Vector)
  if length(grad) > 0
    g!(grad,x)
  end
  integral = integral_function(x)
  global iter
  iter::Int += 1
  if mod(iter,1000) == 0
    global g_res
    g_res[1] = a2*B021 + a3*B031 + a3*B032 - 1
    g_res[2] = a2*B021^2 + a3*(B031 + B032)^2 - 3/2
    g_res[3] = A021*a2 + A031*a3 + A032*a3 - 1/2
    @show x
    @show g_res
    @show integral
  end
  integral
end

function G(A021,A031,A032,a1,a2,a3)
  tmps = zeros(Threads.nthreads())
  xmin = 6
  ymin = 4
  xs = -xmin:0.01:1
  ys = -ymin:0.01:4
  Threads.@threads for x in xs
    for y in ys
      z = x + im*y
      tmps[Threads.threadid()] += abs(1 + z*(a1 + a2 + A021*z*a2 + a3 + A032*z*a3 + (A031*z + A021*A032*z^2)*a3))<1
    end
  end
  tmp = sum(tmps)
  (tmp/(length(xs)*length(ys)))*((xmin+1)*(2ymin))
end

function constraints!(res,x,grad)
  a2, a3, b13, b22, b23, A021, A031, B021, B031, B032 = x
  c01 = 0
  c03 = 1
  c11 = 0
  c13 = 1
  A032 = 1 - A031
  b12 = (b22*(1 - c11 + b13*c11 - b13*c13))/(-1 + b23*c11 - b23*c13)
  b11 = 1 - b12 - b13
  b21 = - b22 - b23
  c12 = (-1 - b21*c11 - b23*c13)/b22
  a1 = 1 - a2 - a3
  c02 = A021
  if length(grad)>0
    grad .= 0.0
    grad[1,1] = B021
    grad[2,1] = B031 + B032
    grad[12,1] = a2
    grad[13,1] = a3
    grad[14,1] = a3
    grad[1,2] = B021^2
    grad[2,2] = (B031 + B032)^2
    grad[12,2] = 2a2*B021
    grad[13,2] = 2a3*(B031+B032)
    grad[14,2] = 2a3*(B031+B032)
    grad[1,3] = A021
    grad[2,3] = A031 + A032
    grad[9,3] = a2
    grad[10,3] = a3
    grad[11,3] = a3
  end
  res[1] = a2*B021 + a3*B031 + a3*B032 - 1
  res[2] = a2*B021^2 + a3*(B031 + B032)^2 - 3/2
  res[3] = A021*a2 + A031*a3 + A032*a3 - 1/2
end

function c_constraints!(res,x,grad)
  a2, a3, b13, b22, b23, A021, A031, B021, B031, B032 = x
  c01 = 0
  c03 = 1
  c11 = 0
  c13 = 1
  A032 = 1 - A031
  b12 = (b22*(1 - c11 + b13*c11 - b13*c13))/(-1 + b23*c11 - b23*c13)
  b11 = 1 - b12 - b13
  b21 = - b22 - b23
  c12 = (-1 - b21*c11 - b23*c13)/b22
  a1 = 1 - a2 - a3
  c02 = A021
  res[1] = c12 - 1
  res[2] = -c12
end

a2 = 1/6; a3 = 2/3; c01 = 0; c11 = 1; c13 = 0; b13 = 0;
b22 = 1; b23 = 0; A021 = 1; A031 = 1/4; A032 = 1/4;
B021 = 0; B031 = 1; B032 = 1/2
#x0 = [a2, a3, c01, c11, c13, b13, b22, b23, A021, A031, A032, B021, B031, B032] + rand(14)
x0 = rand(10)
x0[6] = max.(0,x0[6])
x0[6] = min.(1,x0[6])
opt = Opt(:LN_COBYLA, 10) #:LD_SLSQP, :LN_COBYLA (semi), :GN_ISRES, :LN_AUGLAG_EQ support equality constraints
lower_bounds!(opt, lbs)
upper_bounds!(opt, ubs)
max_objective!(opt, f)
equality_constraint!(opt::Opt,constraints!, [1e-12 for i in 1:3])
inequality_constraint!(opt::Opt,c_constraints!, [1e-12 for i in 1:2])
ftol_abs!(opt,1e-12)
xtol_rel!(opt,1e-12)
maxeval!(opt::Opt, Int(1e6))
maxf,maxx,ret = optimize(opt,x0)
println(maxf)

#maxx_str = "x = [$(maxx[1]),$(maxx[2]),$(maxx[3]),$(maxx[4]),$(maxx[5]),$(maxx[6]),$(maxx[7]),$(maxx[8]),$(maxx[9]),$(maxx[10]),$(maxx[11]),$(maxx[12]),$(maxx[13]),$(maxx[14])]"
#println(maxx_str)

### Result

x = maxx
a2, a3, b13, b22, b23, A021, A031, B021, B031, B032 = x
c01 = 0
c03 = 1
c11 = 0
c13 = 1
A032 = 1 - A031
b12 = (b22*(1 - c11 + b13*c11 - b13*c13))/(-1 + b23*c11 - b23*c13)
b11 = 1 - b12 - b13
b21 = - b22 - b23
c12 = (-1 - b21*c11 - b23*c13)/b22
a1 = 1 - a2 - a3
c02 = A021
G(A021,A031,A032,a1,a2,a3)
g_res = zeros(3); g_res_grad = zeros(14,3); c_res=zeros(4)
constraints!(g_res, x, g_res_grad)
g_res
g_res_grad
c_constraints!(c_res, x, g_res_grad)
c_res

value_str = """
  a1 = $a1;
  a2 = $a2;
  a3 = $a3;
  c01 = $c01;
  c02 = $c02;
  c03 = $c03;
  c11 = $c11;
  c12 = $c12;
  c13 = $c13;
  b11 = $b11;
  b12 = $b12;
  b13 = $b13;
  b21 = $b21;
  b22 = $b22;
  b23 = $b23;
  A021 = $A021;
  A031 = $A031;
  A032 = $A032;
  B021 = $B021;
  B031 = $B031;
  B032 = $B032;
  """
println(value_str)
#=
##################
# Test Integral
G(1,13/56,-1/56,13/4,5/4,-7/2) # 12.063
G(1,1/4,1/4,1/6,1/6,2/3) # 9.0872

###################
# Test constraints

a2 = 1/6; a3 = 2/3; c01 = 0; c11 = 1; c13 = 0; b13 = 0;
b22 = 1; b23 = 0; A021 = 1; A031 = 1/4; A032 = 1/4;
B021 = 0; B031 = 1; B032 = 1/2
x = [a2, a3, c01, c11, c13, b13, b22, b23, A021, A031, A032, B021, B031, B032]


a2, a3, c01, c11, c13, b13, b22, b23, A021, A031, A032, B021, B031, B032 = x
b12 = (b22*(1 - c11 + b13*c11 - b13*c13))/(-1 + b23*c11 - b23*c13)
b11 = 1 - b12 - b13
b21 = - b22 - b23
c12 = (-1 - b21*c11 - b23*c13)/b22
a1 = 1 - a2 - a3
c02 = A021
c03 = A031 + A032

a1 ≈ 1/6
a2 ≈ 1/6
a3 == 2/3
c01 == 0
c02 == 1
c03 == 1/2
c11 == 1
c12 == 0
c13 == 0
A021 == 1
A031 == 1/4
A032 == 1/4
B021 == 0
B031 == 1
B032 == 1/2
b11 == 1
b12 == 0
b13 == 0
b21 == -1
b22 == 1
b23 == 0

g_res = zeros(3); g_res_grad = zeros(14,3); c_res=zeros(4)
constraints!(g_res, x, g_res_grad)
c_constraints!(c_res, x, g_res_grad)
g_res
g_res_grad
c_res
=#
