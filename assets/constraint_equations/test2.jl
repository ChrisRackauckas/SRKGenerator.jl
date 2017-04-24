using JuMP
using GLPKMathProgInterface
using Ipopt

n = 2;
lowRes = -4;
hiRes  = 4;
m = Model()
@defVar(m, lowRes <= α[i=1:n]  <=hiRes)
@defVar(m, lowRes <= β1[i=1:n] <=hiRes)
@defVar(m, lowRes <= β2[i=1:n] <=hiRes)
@defVar(m, lowRes <= β3[i=1:n] <=hiRes)
@defVar(m, lowRes <= β4[i=1:n] <=hiRes)

@defVar(m, lowRes <= A011 <= hiRes)
@defVar(m, lowRes <= A012 <= hiRes)
@defVar(m, lowRes <= A021 <= hiRes)
@defVar(m, lowRes <= A022 <= hiRes)

@defVar(m, lowRes <= A111 <= hiRes)
@defVar(m, lowRes <= A112 <= hiRes)
@defVar(m, lowRes <= A121 <= hiRes)
@defVar(m, lowRes <= A122 <= hiRes)

@defVar(m, lowRes <= B011 <= hiRes)
@defVar(m, lowRes <= B012 <= hiRes)
@defVar(m, lowRes <= B021 <= hiRes)
@defVar(m, lowRes <= B022 <= hiRes)

@defVar(m, lowRes <= B111 <= hiRes)
@defVar(m, lowRes <= B112 <= hiRes)
@defVar(m, lowRes <= B121 <= hiRes)
@defVar(m, lowRes <= B122 <= hiRes)

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
@addConstraints(m,begin
                  sum{α[i], i=1:n}==1
                  sum{β1[i], i=1:n}==1
                  sum{β1[i], i=1:n}==1
                  sum{β2[i], i=1:n}==0
                  sum{β3[i], i=1:n}==0
                  sum{β4[i], i=1:n}==0
                  sum{β4[i], i=1:n}==0
                end)
@addNLConstraints(m,begin
                  B111*β1[1]+B112*β1[1]+B121*β1[2]+B122*β1[2]==0
                  B111*β2[1]+B112*β2[1]+B121*β2[2]+B122*β2[2]==1
                  B111*β3[1]+B112*β3[1]+B121*β3[2]+B122*β3[2]==0
                  B111*β4[1]+B112*β4[1]+B121*β4[2]+B122*β4[2]==0
                  A011*α[1] +A012*α[1] +A021*α[2] +A022*α[2] ==1//2
                  B011*α[1] +B012*α[1] +B021*α[2] +B022*α[2] ==1
                  A111*β1[1]+A112*β1[1]+A121*β1[2]+A122*β1[2]==1
                  A111*β2[1]+A112*β2[1]+A121*β2[2]+A122*β2[2]==0
                  A111*β3[1]+A112*β3[1]+A121*β3[2]+A122*β3[2]==-1
                  A111*β4[1]+A112*β4[1]+A121*β4[2]+A122*β4[2]==0
                end)

@addNLConstraint(m,α[1] * (B011+B012)*(B011+B012)  + α[2] * (B021+B022)*(B021+B022)  == 3//2)
@addNLConstraint(m,β1[1] * (B111+B112)*(B111+B112)  + β1[2] * (B121+B122)*(B121+B122)  == 1)
@addNLConstraint(m,β2[1] * (B111+B112)*(B111+B112)  + β2[2] * (B121+B122)*(B121+B122)  == 0)
@addNLConstraint(m,β3[1] * (B111+B112)*(B111+B112)  + β3[2] * (B121+B122)*(B121+B122)  == -1)
@addNLConstraint(m,β4[1] * (B111+B112)*(B111+B112)  + β4[2] * (B121+B122)*(B121+B122)  == 2)
@addNLConstraint(m,β1[1] * (B111*(B111+B112) + B112*(B121+B122))  + β1[2] * (B121*(B111+B112) + B122*(B121+B122))  == 0)
@addNLConstraint(m,β2[1] * (B111*(B111+B112) + B112*(B121+B122))  + β2[2] * (B121*(B111+B112) + B122*(B121+B122))  == 0)
@addNLConstraint(m,β3[1] * (B111*(B111+B112) + B112*(B121+B122))  + β3[2] * (B121*(B111+B112) + B122*(B121+B122))  == 0)
@addNLConstraint(m,β4[1] * (B111*(B111+B112) + B112*(B121+B122))  + β4[2] * (B121*(B111+B112) + B122*(B121+B122))  == 1)

@addNLConstraint(m,(1/2)*((A111*(B011+B012) + A112*(B021+B022))*β1[1] +
      (A121*(B011+B012) + A122*(B021+B022))β1[2]) +
      (1/3)*((A111*(B011+B012) + A112*(B021+B022))β3[1] +
      (A121*(B011+B012) + A122*(B021+B022))β3[2]) == 0)

print(m)

status = solve(m)

println("Objective value: ", getObjectiveValue(m))
println("β1 = ", getValue(β1))
println("β2 = ", getValue(β2))
println("β3 = ", getValue(β3))
println("β4 = ", getValue(β4))
