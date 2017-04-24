using JuMP
using GLPKMathProgInterface
using Ipopt

n = 2;
m = Model()
@defVar(m, -2 <= α[i=1:n] <= 2)
@defVar(m, -2 <= β1[i=1:n] <= 2)
@defVar(m, -2 <= β2[i=1:n] <= 2)
@defVar(m, -2 <= β3[i=1:n] <= 2)
@defVar(m, -2 <= β4[i=1:n] <= 2)

@defVar(m, -2 <= A011 <= 2)
@defVar(m, -2 <= A012 <= 2)
@defVar(m, -2 <= A021 <= 2)
@defVar(m, -2 <= A022 <= 2)

@defVar(m, -2 <= A111 <= 2)
@defVar(m, -2 <= A112 <= 2)
@defVar(m, -2 <= A121 <= 2)
@defVar(m, -2 <= A122 <= 2)

@defVar(m, -2 <= B011 <= 2)
@defVar(m, -2 <= B012 <= 2)
@defVar(m, -2 <= B021 <= 2)
@defVar(m, -2 <= B022 <= 2)

@defVar(m, -2 <= B111 <= 2)
@defVar(m, -2 <= B112 <= 2)
@defVar(m, -2 <= B121 <= 2)
@defVar(m, -2 <= B122 <= 2)

A0 = [A011 A012
     A021 A022]
A1 = [A111 A112
     A121 A122]
B0 = [B011 B012
      B021 B022]
B1 = [B111 B112
      B121 B122]
e = [1 1]'

@setObjective(m, Max, 0)
@addConstraints(m,sum{α[i], i=1:n}==1)
@addConstraint(m,sum{β1[i], i=1:n}==1)
@addConstraint(m,sum{β2[i], i=1:n}==0)
@addConstraint(m,sum{β3[i], i=1:n}==0)
@addConstraint(m,sum{β4[i], i=1:n}==0)
@addConstraint(m,sum{β4[i], i=1:n}==0)
@addNLConstraint(m,transpose(β1)*B1*e==0)
@addNLConstraint(m,transpose(β2)*B1*e==1)
@addNLConstraint(m,transpose(β3)*B1*e==0)
@addNLConstraint(m,transpose(β4)*B1*e==0)

@addNLConstraint(m,transpose(α)*A0*e==1//2)
@addNLConstraint(m,transpose(α)*B0*e==1)
@addNLConstraint(m,transpose(α)*(B0*e).^2 ==3//2)
@addNLConstraint(m,transpose(β1)*A1*e==1)
@addNLConstraint(m,transpose(β2)*A1*e==0)
@addNLConstraint(m,transpose(β3)*A1*e==-1)
@addNLConstraint(m,transpose(β4)*A1*e==0)

@addNLConstraint(m,transpose(β1)*(B1*e).^2==1)
@addNLConstraint(m,transpose(β3)*(B1*e).^2==-1)
@addNLConstraint(m,transpose(β2)*(B1*e).^2==0)
@addNLConstraint(m,transpose(β4)*(B1*e).^2==2)

@addNLConstraint(m,transpose(β1)*(B1*(B1*e))==0)
@addNLConstraint(m,transpose(β2)*(B1*(B1*e))==0)
@addNLConstraint(m,transpose(β3)*(B1*(B1*e))==0)
@addNLConstraint(m,transpose(β4)*(B1*(B1*e))==1)

@addNLConstraint(m,(1//2)*transpose(β1)*(A1*(B0*e))+(1//3)*transpose(β3)*(A1*(B0*e))==0)

print(m)

status = solve(m)

println("Objective value: ", getObjectiveValue(m))
println("x = ", getValue(x))
println("y = ", getValue(y))
