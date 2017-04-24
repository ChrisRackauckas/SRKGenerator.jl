using JuMP
using GLPKMathProgInterface
using Ipopt

n = 3;
lowRes = -20;
hiRes  = 20;
m = Model(solver=IpoptSolver())
@defVar(m, lowRes <= α[i=1:n]  <=hiRes)
@defVar(m, lowRes <= β1[i=1:n] <=hiRes)
@defVar(m, lowRes <= β2[i=1:n] <=hiRes)
@defVar(m, lowRes <= β3[i=1:n] <=hiRes)
@defVar(m, lowRes <= β4[i=1:n] <=hiRes)

@defVar(m, lowRes <= A011 <= hiRes)
@defVar(m, lowRes <= A012 <= hiRes)
@defVar(m, lowRes <= A013 <= hiRes)
@defVar(m, lowRes <= A021 <= hiRes)
@defVar(m, lowRes <= A022 <= hiRes)
@defVar(m, lowRes <= A023 <= hiRes)
@defVar(m, lowRes <= A031 <= hiRes)
@defVar(m, lowRes <= A032 <= hiRes)
@defVar(m, lowRes <= A033 <= hiRes)

@defVar(m, lowRes <= A111 <= hiRes)
@defVar(m, lowRes <= A112 <= hiRes)
@defVar(m, lowRes <= A113 <= hiRes)
@defVar(m, lowRes <= A121 <= hiRes)
@defVar(m, lowRes <= A122 <= hiRes)
@defVar(m, lowRes <= A123 <= hiRes)
@defVar(m, lowRes <= A131 <= hiRes)
@defVar(m, lowRes <= A132 <= hiRes)
@defVar(m, lowRes <= A133 <= hiRes)

@defVar(m, lowRes <= B011 <= hiRes)
@defVar(m, lowRes <= B012 <= hiRes)
@defVar(m, lowRes <= B013 <= hiRes)
@defVar(m, lowRes <= B021 <= hiRes)
@defVar(m, lowRes <= B022 <= hiRes)
@defVar(m, lowRes <= B023 <= hiRes)
@defVar(m, lowRes <= B031 <= hiRes)
@defVar(m, lowRes <= B032 <= hiRes)
@defVar(m, lowRes <= B033 <= hiRes)

@defVar(m, lowRes <= B111 <= hiRes)
@defVar(m, lowRes <= B112 <= hiRes)
@defVar(m, lowRes <= B113 <= hiRes)
@defVar(m, lowRes <= B121 <= hiRes)
@defVar(m, lowRes <= B122 <= hiRes)
@defVar(m, lowRes <= B123 <= hiRes)
@defVar(m, lowRes <= B131 <= hiRes)
@defVar(m, lowRes <= B132 <= hiRes)
@defVar(m, lowRes <= B133 <= hiRes)

A0 = [A011 A012
     A021 A022]
A1 = [A111 A112
     A121 A122]
B0 = [B011 B012
      B021 B022]
B1 = [B111 B112
      B121 B122]
e = [1 1]'

#@setObjective(m, Max, 0)
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
                  B111*β1[1]+B112*β1[1]+B113*β1[1] + B121*β1[2]+B122*β1[2]+B123*β1[2] + B131*β1[3]+B132*β1[3]+B133*β1[3] ==0
                  B111*β2[1]+B112*β2[1]+B113*β2[1] + B121*β2[2]+B122*β2[2]+B123*β2[2] + B131*β2[3]+B132*β2[3]+B133*β2[3] ==1
                  B111*β3[1]+B112*β3[1]+B113*β3[1] + B121*β3[2]+B122*β3[2]+B123*β3[2] + B131*β3[3]+B132*β3[3]+B133*β3[3] ==0
                  B111*β4[1]+B112*β4[1]+B113*β4[1] + B121*β4[2]+B122*β4[2]+B123*β4[2] + B131*β4[3]+B132*β4[3]+B133*β4[3] ==0
                  A011*α[1]+A012*α[1]+A013*α[1]    + A021*α[2]+A022*α[2]+A023*α[2]    + A031*α[3]+A032*α[3]+A033*α[3] ==1//2
                  B011*α[1]+B012*α[1]+B013*α[1]    + B021*α[2]+B022*α[2]+B023*α[2]    + B031*α[3]+B032*α[3]+B033*α[3] ==1
                  A111*β1[1]+A112*β1[1]+A113*β1[1] + A121*β1[2]+A122*β1[2]+A123*β1[2] + A131*β1[3]+A132*β1[3]+A133*β1[3]==1
                  A111*β2[1]+A112*β2[1]+A113*β2[1] + A121*β2[2]+A122*β2[2]+A123*β2[2] + A131*β2[3]+A132*β2[3]+A133*β2[3]==0
                  A111*β3[1]+A112*β3[1]+A113*β3[1] + A121*β3[2]+A122*β3[2]+A123*β3[2] + A131*β3[3]+A132*β3[3]+A133*β3[3]==-1
                  A111*β4[1]+A112*β4[1]+A113*β4[1] + A121*β4[2]+A122*β4[2]+A123*β4[2] + A131*β4[3]+A132*β4[3]+A133*β4[3]==0
                end)

@addNLConstraint(m,α[1]  * (B011+B012+B013)*(B011+B012+B013)  + α[2]  * (B021+B022+B123)*(B021+B022+B123) + α[3]  * (B031+B032+B133)*(B031+B032+B133)  == 3//2)
@addNLConstraint(m,β1[1] * (B111+B112+B013)*(B111+B112+B013)  + β1[2] * (B121+B122+B123)*(B121+B122+B123) + β1[3] * (B131+B132+B133)*(B131+B132+B133)  == 1)
@addNLConstraint(m,β2[1] * (B111+B112+B013)*(B111+B112+B013)  + β2[2] * (B121+B122+B123)*(B121+B122+B123) + β2[3] * (B131+B132+B133)*(B131+B132+B133)  == 0)
@addNLConstraint(m,β3[1] * (B111+B112+B013)*(B111+B112+B013)  + β3[2] * (B121+B122+B123)*(B121+B122+B123) + β3[3] * (B131+B132+B133)*(B131+B132+B133)  == -1)
@addNLConstraint(m,β4[1] * (B111+B112+B013)*(B111+B112+B013)  + β4[2] * (B121+B122+B123)*(B121+B122+B123) + β4[3] * (B131+B132+B133)*(B131+B132+B133)  == 2)

@addNLConstraint(m,β1[1] * (B111*(B111+B112+B113) + B112*(B121+B122+B123) + B113*(B131+B132+B133))  + β1[2] * (B121*(B111+B112+B113) + B122*(B121+B122+B123) + B123*(B131+B132+B133))  + β1[3] * (B131*(B111+B112+B113) + B132*(B121+B122+B123) + B133*(B131+B132+B133))  == 0)
@addNLConstraint(m,β2[1] * (B111*(B111+B112+B113) + B112*(B121+B122+B123) + B113*(B131+B132+B133))  + β2[2] * (B121*(B111+B112+B113) + B122*(B121+B122+B123) + B123*(B131+B132+B133))  + β2[3] * (B131*(B111+B112+B113) + B132*(B121+B122+B123) + B133*(B131+B132+B133))  == 0)
@addNLConstraint(m,β3[1] * (B111*(B111+B112+B113) + B112*(B121+B122+B123) + B113*(B131+B132+B133))  + β3[2] * (B121*(B111+B112+B113) + B122*(B121+B122+B123) + B123*(B131+B132+B133))  + β3[3] * (B131*(B111+B112+B113) + B132*(B121+B122+B123) + B133*(B131+B132+B133))  == 0)
@addNLConstraint(m,β4[1] * (B111*(B111+B112+B113) + B112*(B121+B122+B123) + B113*(B131+B132+B133))  + β4[2] * (B121*(B111+B112+B113) + B122*(B121+B122+B123) + B123*(B131+B132+B133))  + β4[3] * (B131*(B111+B112+B113) + B132*(B121+B122+B123) + B133*(B131+B132+B133))  == 1)

@addNLConstraint(m,(1/2)*((A111*(B011+B012+B013) + A112*(B021+B022+B023) + A113*(B031+B032+B033))*β1[1] +
      (A121*(B011+B012+B013) + A122*(B021+B022+B023) + A123*(B031+B032+B033))β1[2] +
      (A131*(B011+B012+B013) + A132*(B021+B022+B023) + A133*(B031+B032+B033))β1[3]) +
      (1/3)*((A111*(B011+B012+B013) + A112*(B021+B022+B023) + A113*(B031+B032+B033))*β3[1] +
            (A121*(B011+B012+B013) + A122*(B021+B022+B023) + A123*(B031+B032+B033))β3[2] +
            (A131*(B011+B012+B013) + A132*(B021+B022+B023) + A133*(B031+B032+B033))β3[3]) == 0)

print(m)

status = solve(m)

println("Objective value: ", getObjectiveValue(m))
println("β1 = ", getValue(β1))
println("β2 = ", getValue(β2))
println("β3 = ", getValue(β3))
println("β4 = ", getValue(β4))
