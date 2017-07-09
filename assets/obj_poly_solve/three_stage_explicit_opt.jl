using NLopt, ForwardDiff

x0 = zeros(14)
gcfg = ForwardDiff.GradientConfig(nothing,x0)
g!(grad, x) = ForwardDiff.gradient!(grad, integral_function, x, gcfg)
g_con!(grad, x) = ForwardDiff.gradient!(grad, constraint_cost!, x, gcfg)
lbs  = -5ones(14)
ubs = 20ones(14)

function constraint_cost!(res::Vector, x::Vector)
  a2, a3, c01, c11, c13, b13, b22, b23, A021, A031, A032, B021, B031, B032 = x
  b12 = (b22*(1 - c11 + b13*c11 - b13*c13))/(-1 + b23*c11 - b23*c13)
  b11 = 1 - b12 - b13
  b21 = - b22 - b23
  c12 = (-1 - b21*c11 - b23*c13)/b22
  a1 = 1 - a2 - a3
  c02 = A021
  c03 = A031 + A032
  res[1] = a2*B021 + a3*B031 + a3*B032 - 1
  res[2] = a2*B021^2 + a3*(B031 + B032)^2 - 3/2
end

function integral_function(x)
  a2, a3, c01, c11, c13, b13, b22, b23, A021, A031, A032, B021, B031, B032 = x
  b12 = (b22*(1 - c11 + b13*c11 - b13*c13))/(-1 + b23*c11 - b23*c13)
  b11 = 1 - b12 - b13
  b21 = - b22 - b23
  c12 = (-1 - b21*c11 - b23*c13)/b22
  a1 = 1 - a2 - a3
  c02 = A021
  c03 = A031 + A032
  G(A021,A031,A032,a1,a2,a3)
end

iter = 0
g_res = zeros(2)

function f(x::Vector, grad::Vector)
  if length(grad) > 0
    g!(grad,x)
  end
  integral = integral_function(x)
  global iter
  iter::Int += 1
  if mod(iter,1000) == 0
    global g_res
    constraint_cost!(g_res, x)
    @show x
    @show g_res
    @show integral
  end
  integral
end

function G(A021,A031,A032,a1,a2,a3)
  tmps = zeros(Threads.nthreads())
  xs = -4:0.01:1
  ys = -3:0.01:3
  Threads.@threads for x in xs
    for y in ys
      z = x + im*y
      tmps[Threads.threadid()] += abs(1 + z*(a1 + a2 + A021*z*a2 + a3 + A032*z*a3 + (A031*z + A021*A032*z^2)*a3))<1
    end
  end
  tmp = sum(tmps)
  (tmp/(length(xs)*length(ys)))*(5*6)
end

function constraints!(res,x,grad)
  if length(grad)>0
    g_con!(grad, x)
  end
  constraint_cost!(res, x)
end

a2 = 1/6; a3 = 2/3; c01 = 0; c11 = 1; c13 = 0; b13 = 0;
b22 = 1; b23 = 0; A021 = 1; A031 = 1/4; A032 = 1/4;
B021 = 0; B031 = 1; B032 = 1/2
x0 = [a2, a3, c01, c11, c13, b13, b22, b23, A021, A031, A032, B021, B031, B032]

opt = Opt(:LN_AUGLAG_EQ, 14) #:LD_SLSQP, :LN_COBYLA (semi), :GN_ISRES, :LN_AUGLAG_EQ support equality constraints
lower_bounds!(opt, lbs)
upper_bounds!(opt, ubs)
max_objective!(opt, f)
equality_constraint!(opt::Opt,constraints!, [1e-8 for i in 1:2])
ftol_abs!(opt,1e-8)
xtol_rel!(opt,1e-4)
maxeval!(opt::Opt, Int(1e6))
maxf,maxx,ret = optimize(opt,x0)

### Result

x = [0.208263, 0.666718, -0.0522899, 0.945187, -0.123953, 0.0160279, 0.931976, 0.0433489, 0.903992, 0.117123, 0.0730312, -0.137181, 1.85575, -0.365674]
a2, a3, c01, c11, c13, b13, b22, b23, A021, A031, A032, B021, B031, B032 = x
b12 = (b22*(1 - c11 + b13*c11 - b13*c13))/(-1 + b23*c11 - b23*c13)
b11 = 1 - b12 - b13
b21 = - b22 - b23
c12 = (-1 - b21*c11 - b23*c13)/b22
a1 = 1 - a2 - a3
c02 = A021
c03 = A031 + A032
G(A021,A031,A032,a1,a2,a3)
g_res = [0.0, 0.0]
constraint_cost!(g_res::Vector, x::Vector)
g_res

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

G(A021,A031,A032,a1,a2,a3)
=#
