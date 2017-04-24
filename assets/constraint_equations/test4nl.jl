using JuMP
using GLPKMathProgInterface
using Ipopt
using NLopt

n = 4;
lowRes = -2;
hiRes  = 2;
m = Model(solver=NLoptSolver(algorithm=:LD_SLSQP))
#@setNLObjective(m, Max,0)

@defVar(m, lowRes <= α[i=1:n]  <=hiRes)
@defVar(m, lowRes <= β1[i=1:n] <=hiRes)
@defVar(m, lowRes <= β2[i=1:n] <=hiRes)
@defVar(m, lowRes <= β3[i=1:n] <=hiRes)
@defVar(m, lowRes <= β4[i=1:n] <=hiRes)

@defVar(m, lowRes <= A011 <= hiRes)
@defVar(m, lowRes <= A012 <= hiRes)
@defVar(m, lowRes <= A013 <= hiRes)
@defVar(m, lowRes <= A014 <= hiRes)
@defVar(m, lowRes <= A021 <= hiRes)
@defVar(m, lowRes <= A022 <= hiRes)
@defVar(m, lowRes <= A023 <= hiRes)
@defVar(m, lowRes <= A024 <= hiRes)
@defVar(m, lowRes <= A031 <= hiRes)
@defVar(m, lowRes <= A032 <= hiRes)
@defVar(m, lowRes <= A033 <= hiRes)
@defVar(m, lowRes <= A034 <= hiRes)
@defVar(m, lowRes <= A041 <= hiRes)
@defVar(m, lowRes <= A042 <= hiRes)
@defVar(m, lowRes <= A043 <= hiRes)
@defVar(m, lowRes <= A044 <= hiRes)

@defVar(m, lowRes <= A111 <= hiRes)
@defVar(m, lowRes <= A112 <= hiRes)
@defVar(m, lowRes <= A113 <= hiRes)
@defVar(m, lowRes <= A114 <= hiRes)
@defVar(m, lowRes <= A121 <= hiRes)
@defVar(m, lowRes <= A122 <= hiRes)
@defVar(m, lowRes <= A123 <= hiRes)
@defVar(m, lowRes <= A124 <= hiRes)
@defVar(m, lowRes <= A131 <= hiRes)
@defVar(m, lowRes <= A132 <= hiRes)
@defVar(m, lowRes <= A133 <= hiRes)
@defVar(m, lowRes <= A134 <= hiRes)
@defVar(m, lowRes <= A141 <= hiRes)
@defVar(m, lowRes <= A142 <= hiRes)
@defVar(m, lowRes <= A143 <= hiRes)
@defVar(m, lowRes <= A144 <= hiRes)

@defVar(m, lowRes <= B011 <= hiRes)
@defVar(m, lowRes <= B012 <= hiRes)
@defVar(m, lowRes <= B013 <= hiRes)
@defVar(m, lowRes <= B014 <= hiRes)
@defVar(m, lowRes <= B021 <= hiRes)
@defVar(m, lowRes <= B022 <= hiRes)
@defVar(m, lowRes <= B023 <= hiRes)
@defVar(m, lowRes <= B024 <= hiRes)
@defVar(m, lowRes <= B031 <= hiRes)
@defVar(m, lowRes <= B032 <= hiRes)
@defVar(m, lowRes <= B033 <= hiRes)
@defVar(m, lowRes <= B034 <= hiRes)
@defVar(m, lowRes <= B041 <= hiRes)
@defVar(m, lowRes <= B042 <= hiRes)
@defVar(m, lowRes <= B043 <= hiRes)
@defVar(m, lowRes <= B044 <= hiRes)

@defVar(m, lowRes <= B111 <= hiRes)
@defVar(m, lowRes <= B112 <= hiRes)
@defVar(m, lowRes <= B113 <= hiRes)
@defVar(m, lowRes <= B114 <= hiRes)
@defVar(m, lowRes <= B121 <= hiRes)
@defVar(m, lowRes <= B122 <= hiRes)
@defVar(m, lowRes <= B123 <= hiRes)
@defVar(m, lowRes <= B124 <= hiRes)
@defVar(m, lowRes <= B131 <= hiRes)
@defVar(m, lowRes <= B132 <= hiRes)
@defVar(m, lowRes <= B133 <= hiRes)
@defVar(m, lowRes <= B134 <= hiRes)
@defVar(m, lowRes <= B141 <= hiRes)
@defVar(m, lowRes <= B142 <= hiRes)
@defVar(m, lowRes <= B143 <= hiRes)
@defVar(m, lowRes <= B144 <= hiRes)

#=
A0 = [A011 A012
     A021 A022]
A1 = [A111 A112
     A121 A122]
B0 = [B011 B012
      B021 B022]
B1 = [B111 B112
      B121 B122]
e = [1 1]'
=#

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
                  B111*β1[1]+B112*β1[1]+B113*β1[1]+B114*β1[1] + B121*β1[2]+B122*β1[2]+B123*β1[2]+B124*β1[2] + B131*β1[3]+B132*β1[3]+B133*β1[3]+B134*β1[3]+ B141*β1[4]+B142*β1[4]+B143*β1[4]+B144*β1[4] ==0
                  B111*β2[1]+B112*β2[1]+B113*β2[1]+B114*β2[1] + B121*β2[2]+B122*β2[2]+B123*β2[2]+B124*β2[2] + B131*β2[3]+B132*β2[3]+B133*β2[3]+B134*β2[3]+ B141*β2[4]+B142*β2[4]+B143*β2[4]+B144*β2[4] ==1
                  B111*β3[1]+B112*β3[1]+B113*β3[1]+B114*β3[1] + B121*β3[2]+B122*β3[2]+B123*β3[2]+B124*β3[2] + B131*β3[3]+B132*β3[3]+B133*β3[3]+B134*β3[3]+ B141*β3[4]+B142*β3[4]+B143*β3[4]+B144*β3[4] ==0
                  B111*β4[1]+B112*β4[1]+B113*β4[1]+B114*β4[1] + B121*β4[2]+B122*β4[2]+B123*β4[2]+B124*β4[2] + B131*β4[3]+B132*β4[3]+B133*β4[3]+B134*β4[3]+ B141*β4[4]+B142*β4[4]+B143*β4[4]+B144*β4[4] ==0
                  A011*α[1]+A012*α[1]+A013*α[1]+A014*α[1]     + A021*α[2]+A022*α[2]+A023*α[2]+A024*α[2]    + A031*α[3]+A032*α[3]+A033*α[3]+A034*α[3] + A041*α[4]+A042*α[4]+A043*α[4]+A044*α[4] ==1//2
                  B011*α[1]+B012*α[1]+B013*α[1]+B014*α[1]     + B021*α[2]+B022*α[2]+B023*α[2]+B024*α[2]    + B031*α[3]+B032*α[3]+B033*α[3]+B034*α[3] + B041*α[4]+B042*α[4]+B043*α[4]+B044*α[4] ==1
                  A111*β1[1]+A112*β1[1]+A113*β1[1]+A114*β1[1] + A121*β1[2]+A122*β1[2]+A123*β1[2]+A124*β1[2] + A131*β1[3]+A132*β1[3]+A133*β1[3]+A134*β1[3] + A141*β1[4]+A142*β1[4]+A143*β1[4]+A144*β1[4]==1
                  A111*β2[1]+A112*β2[1]+A113*β2[1]+A114*β2[1] + A121*β2[2]+A122*β2[2]+A123*β2[2]+A124*β2[2] + A131*β2[3]+A132*β2[3]+A133*β2[3]+A134*β2[3] + A141*β2[4]+A142*β2[4]+A143*β2[4]+A144*β2[4]==0
                  A111*β3[1]+A112*β3[1]+A113*β3[1]+A114*β3[1] + A121*β3[2]+A122*β3[2]+A123*β3[2]+A124*β3[2] + A131*β3[3]+A132*β3[3]+A133*β3[3]+A134*β3[3] + A141*β3[4]+A142*β3[4]+A143*β3[4]+A144*β3[4]==-1
                  A111*β4[1]+A112*β4[1]+A113*β4[1]+A114*β4[1] + A121*β4[2]+A122*β4[2]+A123*β4[2]+A124*β4[2] + A131*β4[3]+A132*β4[3]+A133*β4[3]+A134*β4[3] + A141*β4[4]+A142*β4[4]+A143*β4[4]+A144*β4[4]==0
                end)

@addNLConstraint(m,α[1]  * (B011+B012+B013+B014)*(B011+B012+B013+B014)  + α[2]  * (B021+B022+B123+B124)*(B021+B022+B123+B124) + α[4]  * (B041+B042+B143+B144)*(B041+B042+B143+B144)  == 3//2)
@addNLConstraint(m,β1[1] * (B111+B112+B013+B014)*(B111+B112+B013+B014)  + β1[2] * (B121+B122+B123+B124)*(B121+B122+B123+B124) + β1[4] * (B141+B142+B143+B144)*(B141+B142+B143+B144)  == 1)
@addNLConstraint(m,β2[1] * (B111+B112+B013+B014)*(B111+B112+B013+B014)  + β2[2] * (B121+B122+B123+B124)*(B121+B122+B123+B124) + β2[4] * (B141+B142+B143+B144)*(B141+B142+B143+B144)  == 0)
@addNLConstraint(m,β3[1] * (B111+B112+B013+B014)*(B111+B112+B013+B014)  + β3[2] * (B121+B122+B123+B124)*(B121+B122+B123+B124) + β3[4] * (B141+B142+B143+B144)*(B141+B142+B143+B144)  == -1)
@addNLConstraint(m,β4[1] * (B111+B112+B013+B014)*(B111+B112+B013+B014)  + β4[2] * (B121+B122+B123+B124)*(B121+B122+B123+B124) + β4[4] * (B141+B142+B143+B144)*(B141+B142+B143+B144)  == 2)

@addNLConstraint(m,β1[1] * (B111*(B111+B112+B113+B114) + B112*(B121+B122+B123+B124) + B113*(B131+B132+B133+B134) + B114*(B141+B142+B143+B144))  + β1[2] * (B121*(B111+B112+B113+B114) + B122*(B121+B122+B123+B124) + B123*(B131+B132+B133+B134) + B124*(B141+B142+B143+B144) )  + β1[3] * (B131*(B111+B112+B113+B114) + B132*(B121+B122+B123+B124) + B133*(B131+B132+B133+B134) + B134*(B141+B142+B143+B144)) + β1[4] * (B141*(B111+B112+B113+B114) + B142*(B121+B122+B123+B124) + B143*(B131+B132+B133+B134) + B144*(B141+B142+B143+B144)) == 0)
@addNLConstraint(m,β2[1] * (B111*(B111+B112+B113+B114) + B112*(B121+B122+B123+B124) + B113*(B131+B132+B133+B134) + B114*(B141+B142+B143+B144))  + β2[2] * (B121*(B111+B112+B113+B114) + B122*(B121+B122+B123+B124) + B123*(B131+B132+B133+B134) + B124*(B141+B142+B143+B144) )  + β2[3] * (B131*(B111+B112+B113+B114) + B132*(B121+B122+B123+B124) + B133*(B131+B132+B133+B134) + B134*(B141+B142+B143+B144)) + β2[4] * (B141*(B111+B112+B113+B114) + B142*(B121+B122+B123+B124) + B143*(B131+B132+B133+B134) + B144*(B141+B142+B143+B144)) == 0)
@addNLConstraint(m,β3[1] * (B111*(B111+B112+B113+B114) + B112*(B121+B122+B123+B124) + B113*(B131+B132+B133+B134) + B114*(B141+B142+B143+B144))  + β3[2] * (B121*(B111+B112+B113+B114) + B122*(B121+B122+B123+B124) + B123*(B131+B132+B133+B134) + B124*(B141+B142+B143+B144) )  + β3[3] * (B131*(B111+B112+B113+B114) + B132*(B121+B122+B123+B124) + B133*(B131+B132+B133+B134) + B134*(B141+B142+B143+B144)) + β3[4] * (B141*(B111+B112+B113+B114) + B142*(B121+B122+B123+B124) + B143*(B131+B132+B133+B134) + B144*(B141+B142+B143+B144)) == 0)
@addNLConstraint(m,β4[1] * (B111*(B111+B112+B113+B114) + B112*(B121+B122+B123+B124) + B113*(B131+B132+B133+B134) + B114*(B141+B142+B143+B144))  + β4[2] * (B121*(B111+B112+B113+B114) + B122*(B121+B122+B123+B124) + B123*(B131+B132+B133+B134) + B124*(B141+B142+B143+B144) )  + β4[3] * (B131*(B111+B112+B113+B114) + B132*(B121+B122+B123+B124) + B133*(B131+B132+B133+B134) + B134*(B141+B142+B143+B144)) + β4[4] * (B141*(B111+B112+B113+B114) + B142*(B121+B122+B123+B124) + B143*(B131+B132+B133+B134) + B144*(B141+B142+B143+B144)) == 1)

@addNLConstraint(m,(1/2)*((A111*(B011+B012+B013+B014) + A112*(B021+B022+B023+B024) + A113*(B031+B032+B033+B034) + A114*(B041+B042+B043+B044))*β1[1] +
                          (A121*(B011+B012+B013+B014) + A122*(B021+B022+B023+B024) + A123*(B031+B032+B033+B034) + A124*(B041+B042+B043+B044))*β1[2] +
                          (A131*(B011+B012+B013+B014) + A132*(B021+B022+B023+B024) + A133*(B031+B032+B033+B034) + A134*(B041+B042+B043+B044))*β1[3] +
                          (A141*(B011+B012+B013+B014) + A142*(B021+B022+B023+B024) + A143*(B031+B032+B033+B034) + A144*(B041+B042+B043+B044))*β1[4]) +
                   (1/3)*((A111*(B011+B012+B013+B014) + A112*(B021+B022+B023+B024) + A113*(B031+B032+B033+B034) + A114*(B041+B042+B043+B044))*β3[1] +
                          (A121*(B011+B012+B013+B014) + A122*(B021+B022+B023+B024) + A123*(B031+B032+B033+B034) + A124*(B041+B042+B043+B044))*β3[2] +
                          (A131*(B011+B012+B013+B014) + A132*(B021+B022+B023+B024) + A133*(B031+B032+B033+B034) + A134*(B041+B042+B043+B044))*β3[3] +
                          (A141*(B011+B012+B013+B014) + A142*(B021+B022+B023+B024) + A143*(B031+B032+B033+B034) + A144*(B041+B042+B043+B044))*β3[4]) == 0)

print(m)

status = solve(m)

println("Objective value: ", getObjectiveValue(m))
println("α1 = ", getValue(α))
println("β1 = ", getValue(β1))
println("β2 = ", getValue(β2))
println("β3 = ", getValue(β3))
println("β4 = ", getValue(β4))
println("A111 = ", getValue(A111))
println("A121 = ", getValue(A121))
