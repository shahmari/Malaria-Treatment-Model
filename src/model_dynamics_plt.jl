using Plots, DifferentialEquations, LaTeXStrings, Plots.Measures

Plots.default(titlefontsize=16, tickfontsize=14, labelfontsize=16, legendfontsize=10, rightmargin=5mm,
    fontfamily="Sans-serif", frame=:box, label=nothing)

include("basic_reproduction_number.jl")

using .BasicReproductionNumber

FigureDir = "../fig/" # Directory for storing the figures
cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file
println("Directories and modules have been set...\n")

# Define a function that returns a plausible parameter set
function default_parameters()
    return Parameters(
        0.2,    # Biting rate
        0.5,    # Human-to-mosquito transmission
        0.5,    # Mosquito-to-human transmission
        20.0,   # Mosquito-to-human ratio
        0.01,    # Human recovery rate
        0.12,   # Mosquito death rate
        0.1,   # Treatment waning rate
        0.1, # Treatment encounter rate
        0.2,  # Untreated E1 -> E2 rate
        0.2,  # Untreated E2 -> I rate
        0.1, # Treated E1 -> E2 rate
        0.1   # Treated E2 -> I rate
    )
end

u0 = [0.99, 0.01, 0.49, 0.01, 0.0, 0.0, 0.49, 0.01, 0.0, 0.0]
p = default_parameters()

tspan = (0.0, 100.0)
prob = ODEProblem(malaria_ode!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
symbolic_EE = compute_endemic_equilibrium(p)

# Plot the solution
plot(sol, vars=(0, 1), label=L"S_H", lw=2, c=:blue)
plot!(sol, vars=(0, 2), label=L"I_H", lw=2, c=:red)
plot!(sol, vars=(0, 3), label=L"S_M", c=:black)
plot!(sol, vars=(0, 4), label=L"E_{1,M}", c=:green)
plot!(sol, vars=(0, 5), label=L"E_{2,M}", c=:orange)
plot!(sol, vars=(0, 6), label=L"I_M", c=:purple)
plot!(sol, vars=(0, 7), label=L"S_T", ls=:dash, c=:black)
plot!(sol, vars=(0, 8), label=L"E_{1,T}", ls=:dash, c=:green)
plot!(sol, vars=(0, 9), label=L"E_{2,T}", ls=:dash, c=:orange)
plot!(sol, vars=(0, 10), label=L"I_T", ls=:dash, c=:purple)
hline!([symbolic_EE[1]], label=L"S_H^*", c=:blue, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[2]], label=L"I_H^*", c=:red, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[3]], label=L"S_M^*", c=:black, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[4]], label=L"E_{1,M}^*", c=:green, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[5]], label=L"E_{2,M}^*", c=:orange, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[6]], label=L"I_M^*", c=:purple, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[7]], label=L"S_T^*", c=:black, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[8]], label=L"E_{1,T}^*", c=:green, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[9]], label=L"E_{2,T}^*", c=:orange, ls=:dashdot, lw=0.5)
hline!([symbolic_EE[10]], label=L"I_T^*", c=:purple, ls=:dashdot, lw=0.5)
plot!(xlabel="Time (days)", ylabel="Frequency", title="Model dynamics", legend=:outerright, size=(700, 570))

# Save the figure
savefig(FigureDir * "model_dynamics.png")
savefig(FigureDir * "model_dynamics.pdf")
println("Model dynamics have been plotted and saved...\n")