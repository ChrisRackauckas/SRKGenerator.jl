using StochasticDiffEq
using Plots; pyplot()
using LaTeXStrings
f(u,p,t) = -1000.0u*(1-u)*(2-u)
g(u,p,t) = 10.0
prob = SDEProblem(f,g,2.0,(0.0,5.0))

sol = solve(prob,SOSRI(),abstol=1e-2,reltol=1e-2,saveat=0.025,seed=3)
plot(sol,legend=false,yaxis = L"X_t",
     title = "Pathwise Stiffness Example",
     xtickfont = font(16, "Latin Modern Roman"),
     titlefont = font(20, "Latin Modern Roman"),
     ytickfont = font(16, "Latin Modern Roman"),
     guidefont = font(18, "Latin Modern Roman"))

savefig("pathwise_stiffness.png")
savefig("pathwise_stiffness.pdf")
