# scripts/generate_figures.jl

"""
This script serves as the main driver for the malaria transmission modeling project.
It uses the MalariaTransmissionModel module to perform simulations and calculations,
and generates all the figures required by the project outline.
"""

# --- 1. Project Setup and Dependencies ---

# The Pkg module is used to manage Julia's environments.
cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file
using Pkg

# Activate the project environment located in the parent directory ("..").
# This ensures that the script uses the exact package versions defined in
# the project's Project.toml and Manifest.toml files for reproducibility.
Pkg.activate("..")

# Load the necessary packages.
using Plots
using LaTeXStrings # For using LaTeX formatting in plot labels
using DifferentialEquations
using Printf

# Load the core model logic from the `src` directory.
include("../src/MalariaTransmissionModel.jl")
using .MalariaTransmissionModel


# --- 2. Global Configuration & Path Settings ---

# Define a constant for the output directory.
const PLOTS_DIR = "../plots"

# Create the output directory if it doesn't already exist.
if !isdir(PLOTS_DIR)
    mkdir(PLOTS_DIR)
    println("Created directory: '$PLOTS_DIR'")
end

println("Plots will be saved in the '$PLOTS_DIR' directory.")

# Set a consistent, publication-quality theme for all plots.
default(
    fontfamily="Sans-serif",
    linewidth=2,
    size=(800, 600),
    titlefontsize=14,
    guidefontsize=12,
    tickfontsize=10,
    legendfontsize=10,
    framestyle=:box,
    dpi=300
)

println("\nProject environment and settings configured successfully.")
println("Ready to generate figures.")


# --- 3. Helper Functions ---

"""
    setup_scenario_parameters(R0_baseline, a_scale)

A helper function to create a ModelParameters object for a given transmission scenario.
The formulas are based on the project outline.[5]
"""
function setup_scenario_parameters(R0_baseline::Float64, a_scale::Float64)
    # These base parameters are fixed across scenarios as per the outline [5]
    b = 0.5
    c = 0.5
    r = 0.01
    g = 0.05
    L = 2
    theta = 2.0
    s_M = 1.0 / 14.0
    s_T = s_M / theta
    h1 = 0.0 # No waning for this figure
    h2 = 0.0

    # Scenario-dependent parameters [5]
    a = sqrt(R0_baseline) * 0.01 * a_scale
    m = 20.0 / (a_scale^2)

    return ModelParameters(r, a, m, g, L, s_M, s_T, b, c, h1, h2, theta)
end


# --- 4. Figure Generation ---

"""
    generate_figure3()

Generates and saves the plot for Figure 3, showing the impact of introducing
treated nets. The simulation runs in two stages: pre- and post-intervention.
"""
function generate_figure3()
    println("Generating Figure 3: Model Dynamics...")

    # 1. Set up the parameters for a representative scenario.
    params = setup_scenario_parameters(5.0, 2.0)

    # --- STAGE 1: SIMULATE THE RISE OF THE EPIDEMIC (NO NETS) ---

    # 2. Define the initial conditions: a mostly susceptible population with a
    u0_stage1 = [0.9, 0.1, 0.8, 0.05, 0.05, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]
    tspan_stage1 = (0.0, 500.0)  # Simulate for 200 days
    p_stage1 = (params, 0.0, 0.0) # No intervention (p=0, ε=0)
    prob_stage1 = ODEProblem(malaria_ode!, u0_stage1, tspan_stage1, p_stage1)
    sol_stage1 = solve(prob_stage1, Tsit5(), reltol=1e-6, abstol=1e-6)
    ee_sol_stage1 = compute_endemic_equilibrium(params, 0.0, 0.0)
    println("Stage 1 simulation complete.")

    u0_stage2 = hcat(sol_stage1.u...)'[end, :] # Use the final state from stage 1 as initial condition for stage 2
    # --- STAGE 2: INTRODUCE INTERVENTION (NETS) AND SIMULATE DYNAMICS ---
    # 3. Define the intervention parameters and simulate the dynamics post-intervention.
    p_intervention = 0.5      # 50% coverage with insecticide-treated nets
    epsilon_intervention = 0.5 # 50% efficacy of the nets
    t_intervention = 500.0    # Time of intervention (end of stage 1)
    tspan_stage2 = (t_intervention, t_intervention + 750.0) # Simulate for another 200 days
    p_stage2 = (params, p_intervention, epsilon_intervention)
    prob_stage2 = ODEProblem(malaria_ode!, u0_stage2, tspan_stage2, p_stage2)
    sol_stage2 = solve(prob_stage2, Tsit5(), reltol=1e-6, abstol=1e-6)
    println("Stage 2 simulation complete.")

    # --- PLOTTING ---
    # Plot the results from Stage 1
    final_plot = begin
        # --- 3. Figure Generation ---
        # 3a. Plot the epidemic curve from Stage 1 to visualize the outbreak dynamics.
        plot(
            sol_stage1.t,
            hcat(sol_stage1.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            title="Figure 3: Model Dynamics",
            xlabel="Time (days)",
            ylabel="Proportion of Population",
            labels=[L"S_H" L"I_H" L"S_M" L"I_M"],
            color=[:blue :red :green :orange],
            legend=:right
        )
        # for (i, y) in enumerate([ee_sol_stage1.S_H, ee_sol_stage1.I_H, ee_sol_stage1.S_M, ee_sol_stage1.I_M])
        #     plot!([-10, 500], [y, y], lw=0.5, color=([:blue, :red, :green, :orange][i]), label=([L"S_H^*", L"I_H^*", L"S_M^*", L"I_M^*"][i]))
        # end

        plot!(
            sol_stage2.t,
            hcat(sol_stage2.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            color=[:blue :red :green :orange],
            label=["" "" "" ""]
        )

        vline!([t_intervention], lw=2, ls=:dash, color=:black, label="Intervention Start")
    end

    # Save the completed figure
    savefig(final_plot, joinpath(PLOTS_DIR, "figure3_intervention_dynamics.png"))
    println("... Figure 3 saved successfully.")
end


# --- 5. Main Execution Block ---

"""
The main function to run our script.
"""
function main()
    println("\n--- Starting Figure Generation ---")
    generate_figure3()
    # Calls for other figures will be added here.
    println("\n--- Figure Generation Complete ---")
end

# This ensures that main() is called only when the script is executed directly.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end