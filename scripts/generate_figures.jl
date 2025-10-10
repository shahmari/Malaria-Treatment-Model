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
using Plots, Plots.Measures
using LaTeXStrings # For using LaTeX formatting in plot labels
using DifferentialEquations

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
    label=nothing,
    leftmargin=5mm,
    rightmargin=5mm,
    dpi=300
)

println("\nProject environment and settings configured successfully.")
println("Ready to generate figures.")

"""
    setup_scenario_parameters(R0_baseline, a_scale)

A helper function to create a ModelParameters object for a given transmission scenario.
The formulas are based on the project outline.[5]
"""
function setup_scenario_parameters(R0_baseline::Float64, a_scale::Float64)
    # These base parameters are fixed across scenarios as per the outline [1]
    b = 0.5
    c = 0.5
    r = 0.01
    g = 0.05
    L = 2
    theta = 2.0
    s_M = 1.0 / 14.0
    s_T = s_M / theta
    h1 = 0.0
    h2 = 0.0
    a = sqrt(R0_baseline) * 0.01 * a_scale
    m = 20.0 / (a_scale^2)
    p = 0.5
    epsilon = 0.5
    latent_treatment_switch = 0.0

    return ModelParameters(r, a, m, g, L, s_M, s_T, b, c, h1, h2, theta, p, epsilon, latent_treatment_switch)
end


# --- 4. Figure Generation ---

"""
    generate_figure3()

Generates and saves the plot for Figure 3, showing the impact of introducing
treated nets. The simulation runs in two stages: pre- and post-intervention.
"""
function generate_figure3()
    println("Generating Figure 3: Model Dynamics...")

    # Scenario Parameters: R0_baseline ∈ [1.2, 5, 10] && a_scale ∈ [1, 2]
    # Scenario [1, 1]: R0_baseline = 1.2, a_scale = 0.01
    plot_1_1 = begin
        params = setup_scenario_parameters(1.2, 1.0)
        # --- STAGE 1: SIMULATE THE RISE OF THE EPIDEMIC (NO NETS) ---
        u0_stage1 = [0.9, 0.1, 0.8, 0.05, 0.05, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]
        tspan_stage1 = (0.0, 500.0)  # Simulate for 200 days
        params.p = 0.0      # No nets in stage 1
        params.epsilon = 0.0 # No nets in stage 1
        prob_stage1 = ODEProblem(malaria_ode!, u0_stage1, tspan_stage1, params)
        sol_stage1 = solve(prob_stage1, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- STAGE 2: INTRODUCE INTERVENTION (NETS) AND SIMULATE DYNAMICS ---
        u0_stage2 = hcat(sol_stage1.u...)'[end, :] # Use the final state from stage 1 as initial condition for stage 2
        params.p = 0.5      # 50% coverage with insecticide-treated nets
        params.epsilon = 0.5 # 50% efficacy of the nets
        t_intervention = 500.0    # Time of intervention (end of stage 1)
        tspan_stage2 = (t_intervention, t_intervention + 750.0)
        prob_stage2 = ODEProblem(malaria_ode!, u0_stage2, tspan_stage2, params)
        sol_stage2 = solve(prob_stage2, Tsit5(), reltol=1e-6, abstol=1e-6)

        # --- PLOTTING ---
        plot(
            sol_stage1.t,
            hcat(sol_stage1.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            title="Low Prevalence, Low Biting Rate",
            xlabel="Time (days)",
            ylabel="Proportion of Population",
            color=[:blue :red :green :orange],
            legend=:right
        )

        plot!(
            sol_stage2.t,
            hcat(sol_stage2.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            color=[:blue :red :green :orange],
            label=["" "" "" ""]
        )

        vline!([t_intervention], lw=2, ls=:dash, color=:black)
    end
    println("Panel (1,1) done... ")

    # Scenario [1, 2]: R0_baseline = 1.2, a_scale = 2.0
    plot_1_2 = begin
        params = setup_scenario_parameters(1.2, 2.0)
        # --- STAGE 1: SIMULATE THE RISE OF THE EPIDEMIC (NO NETS) ---
        u0_stage1 = [0.9, 0.1, 0.8, 0.05, 0.05, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]
        tspan_stage1 = (0.0, 500.0)  # Simulate for 200 days
        params.p = 0.0      # No nets in stage 1
        params.epsilon = 0.0 # No nets in stage 1
        prob_stage1 = ODEProblem(malaria_ode!, u0_stage1, tspan_stage1, params)
        sol_stage1 = solve(prob_stage1, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- STAGE 2: INTRODUCE INTERVENTION (NETS) AND SIMULATE DYNAMICS ---
        u0_stage2 = hcat(sol_stage1.u...)'[end, :] # Use the final state from stage 1 as initial condition for stage 2
        params.p = 0.5      # 50% coverage with insecticide-treated nets
        params.epsilon = 0.5 # 50% efficacy of the nets
        t_intervention = 500.0    # Time of intervention (end of stage 1)
        tspan_stage2 = (t_intervention, t_intervention + 750.0)
        prob_stage2 = ODEProblem(malaria_ode!, u0_stage2, tspan_stage2, params)
        sol_stage2 = solve(prob_stage2, Tsit5(), reltol=1e-6, abstol=1e-6) # Increased accuracy
        # --- PLOTTING ---
        plot(
            sol_stage1.t,
            hcat(sol_stage1.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            title="Low Prevalence, High Biting Rate",
            xlabel="Time (days)",
            ylabel="Proportion of Population",
            labels=[L"S_H" L"I_H" L"S_M" L"I_M"],
            color=[:blue :red :green :orange],
            legend=:right
        )
        plot!(
            sol_stage2.t,
            hcat(sol_stage2.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            color=[:blue :red :green :orange],
            label=["" "" "" ""]
        )
        vline!([t_intervention], lw=2, ls=:dash, color=:black)
    end
    println("Panel (1,2) done... ")

    # Scenario [2, 1]: R0_baseline = 5.0, a_scale = 1.0
    plot_2_1 = begin
        params = setup_scenario_parameters(5.0, 1.0)
        # --- STAGE 1: SIMULATE THE RISE OF THE EPIDEMIC (NO NETS) ---
        u0_stage1 = [0.9, 0.1, 0.8, 0.05, 0.05, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]
        tspan_stage1 = (0.0, 500.0)  # Simulate for 200 days
        params.p = 0.0      # No nets in stage 1
        params.epsilon = 0.0 # No nets in stage 1
        prob_stage1 = ODEProblem(malaria_ode!, u0_stage1, tspan_stage1, params)
        sol_stage1 = solve(prob_stage1, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- STAGE 2: INTRODUCE INTERVENTION (NETS) AND SIMULATE DYNAMICS ---
        u0_stage2 = hcat(sol_stage1.u...)'[end, :] # Use the final state from stage 1 as initial condition for stage 2
        params.p = 0.5      # 50% coverage with insecticide-treated nets
        params.epsilon = 0.5 # 50% efficacy of the nets
        t_intervention = 500.0    # Time of intervention (end of stage 1)
        tspan_stage2 = (t_intervention, t_intervention + 750.0)
        prob_stage2 = ODEProblem(malaria_ode!, u0_stage2, tspan_stage2, params)
        sol_stage2 = solve(prob_stage2, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- PLOTTING ---
        plot(
            sol_stage1.t,
            hcat(sol_stage1.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            title="Moderate Prevalence, Low Biting Rate",
            xlabel="Time (days)",
            ylabel="Proportion of Population",
            color=[:blue :red :green :orange],
            legend=:right
        )
        plot!(
            sol_stage2.t,
            hcat(sol_stage2.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            color=[:blue :red :green :orange],
            label=["" "" "" ""]
        )
        vline!([t_intervention], lw=2, ls=:dash, color=:black)
    end
    println("Panel (2,1) done... ")

    # Scenario [2, 2]: R0_baseline = 5.0, a_scale = 2.0
    plot_2_2 = begin
        params = setup_scenario_parameters(5.0, 2.0)
        # --- STAGE 1: SIMULATE THE RISE OF THE EPIDEMIC (NO NETS) ---
        u0_stage1 = [0.9, 0.1, 0.8, 0.05, 0.05, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]
        tspan_stage1 = (0.0, 500.0)  # Simulate for 200 days
        params.p = 0.0      # No nets in stage 1
        params.epsilon = 0.0 # No nets in stage 1
        prob_stage1 = ODEProblem(malaria_ode!, u0_stage1, tspan_stage1, params)
        sol_stage1 = solve(prob_stage1, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- STAGE 2: INTRODUCE INTERVENTION (NETS) AND SIMULATE DYNAMICS ---
        u0_stage2 = hcat(sol_stage1.u...)'[end, :] # Use the final state from stage 1 as initial condition for stage 2
        params.p = 0.5      # 50% coverage with insecticide-treated nets
        params.epsilon = 0.5 # 50% efficacy of the nets
        t_intervention = 500.0    # Time of intervention (end of stage 1)
        tspan_stage2 = (t_intervention, t_intervention + 750.0)
        prob_stage2 = ODEProblem(malaria_ode!, u0_stage2, tspan_stage2, params)
        sol_stage2 = solve(prob_stage2, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- PLOTTING ---
        plot(
            sol_stage1.t,
            hcat(sol_stage1.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            title="Moderate Prevalence, High Biting Rate",
            xlabel="Time (days)",
            ylabel="Proportion of Population",
            color=[:blue :red :green :orange],
            legend=:right
        )
        plot!(
            sol_stage2.t,
            hcat(sol_stage2.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            color=[:blue :red :green :orange],
            label=["" "" "" ""]
        )
        vline!([t_intervention], lw=2, ls=:dash, color=:black)
    end
    println("Panel (2,2) done... ")

    # Scenario [3, 1]: R0_baseline = 10.0, a_scale = 1.0
    plot_3_1 = begin
        params = setup_scenario_parameters(10.0, 1.0)
        # --- STAGE 1: SIMULATE THE RISE OF THE EPIDEMIC (NO NETS) ---
        u0_stage1 = [0.9, 0.1, 0.8, 0.05, 0.05, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]
        tspan_stage1 = (0.0, 500.0)  # Simulate for 200 days
        params.p = 0.0      # No nets in stage 1
        params.epsilon = 0.0 # No nets in stage 1
        prob_stage1 = ODEProblem(malaria_ode!, u0_stage1, tspan_stage1, params)
        sol_stage1 = solve(prob_stage1, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- STAGE 2: INTRODUCE INTERVENTION (NETS) AND SIMULATE DYNAMICS ---
        u0_stage2 = hcat(sol_stage1.u...)'[end, :] # Use the final state from stage 1 as initial condition for stage 2
        params.p = 0.5      # 50% coverage with insecticide-treated nets
        params.epsilon = 0.5 # 50% efficacy of the nets
        t_intervention = 500.0    # Time of intervention (end of stage 1)
        tspan_stage2 = (t_intervention, t_intervention + 750.0)
        prob_stage2 = ODEProblem(malaria_ode!, u0_stage2, tspan_stage2, params)
        sol_stage2 = solve(prob_stage2, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- PLOTTING ---
        plot(
            sol_stage1.t,
            hcat(sol_stage1.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            title="High Prevalence, Low Biting Rate",
            xlabel="Time (days)",
            ylabel="Proportion of Population",
            color=[:blue :red :green :orange],
            legend=:right
        )
        plot!(
            sol_stage2.t,
            hcat(sol_stage2.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            color=[:blue :red :green :orange],
            label=["" "" "" ""]
        )
        vline!([t_intervention], lw=2, ls=:dash, color=:black)
    end
    println("Panel (3,1) done... ")

    # Scenario [3, 2]: R0_baseline = 10.0, a_scale = 2.0
    plot_3_2 = begin
        params = setup_scenario_parameters(10.0, 2.0)
        # --- STAGE 1: SIMULATE THE RISE OF THE EPIDEMIC (NO NETS) ---
        u0_stage1 = [0.9, 0.1, 0.8, 0.05, 0.05, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]
        tspan_stage1 = (0.0, 500.0)  # Simulate for 200 days
        params.p = 0.0      # No nets in stage 1
        params.epsilon = 0.0 # No nets in stage 1
        prob_stage1 = ODEProblem(malaria_ode!, u0_stage1, tspan_stage1, params)
        sol_stage1 = solve(prob_stage1, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- STAGE 2: INTRODUCE INTERVENTION (NETS) AND SIMULATE DYNAMICS ---
        u0_stage2 = hcat(sol_stage1.u...)'[end, :] # Use the final state from stage 1 as initial condition for stage 2
        params.p = 0.5      # 50% coverage with insecticide-treated nets
        params.epsilon = 0.5 # 50% efficacy of the nets
        t_intervention = 500.0    # Time of intervention (end of stage 1)
        tspan_stage2 = (t_intervention, t_intervention + 750.0)
        prob_stage2 = ODEProblem(malaria_ode!, u0_stage2, tspan_stage2, params)
        sol_stage2 = solve(prob_stage2, Tsit5(), reltol=1e-6, abstol=1e-6)
        # --- PLOTTING ---
        plot(
            sol_stage1.t,
            hcat(sol_stage1.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            title="High Prevalence, High Biting Rate",
            xlabel="Time (days)",
            ylabel="Proportion of Population",
            color=[:blue :red :green :orange],
            legend=:right
        )
        plot!(
            sol_stage2.t,
            hcat(sol_stage2.u...)'[:, [1, 2, 3, 6]], # Plot S_H, I_H, S_M, I_M
            color=[:blue :red :green :orange],
            label=["" "" "" ""]
        )
        vline!([t_intervention], lw=2, ls=:dash, color=:black)
    end
    println("Panel (3,2) done... ")

    # Combine all subplots into a 3x2 grid layout
    final_plot = plot(
        plot_1_1, plot_1_2,
        plot_2_1, plot_2_2,
        plot_3_1, plot_3_2,
        layout = (3, 2),
        size = (1800, 1600)
    )

    # Save plots
    savefig(plot_1_1, joinpath(PLOTS_DIR, "figure3_panel_1_1.png"))
    savefig(plot_1_1, joinpath(PLOTS_DIR, "figure3_panel_1_1.pdf"))
    savefig(plot_1_2, joinpath(PLOTS_DIR, "figure3_panel_1_2.png"))
    savefig(plot_1_2, joinpath(PLOTS_DIR, "figure3_panel_1_2.pdf"))
    savefig(plot_2_1, joinpath(PLOTS_DIR, "figure3_panel_2_1.png"))
    savefig(plot_2_1, joinpath(PLOTS_DIR, "figure3_panel_2_1.pdf"))
    savefig(plot_2_2, joinpath(PLOTS_DIR, "figure3_panel_2_2.png"))
    savefig(plot_2_2, joinpath(PLOTS_DIR, "figure3_panel_2_2.pdf"))
    savefig(plot_3_1, joinpath(PLOTS_DIR, "figure3_panel_3_1.png"))
    savefig(plot_3_1, joinpath(PLOTS_DIR, "figure3_panel_3_1.pdf"))
    savefig(plot_3_2, joinpath(PLOTS_DIR, "figure3_panel_3_2.png"))
    savefig(plot_3_2, joinpath(PLOTS_DIR, "figure3_panel_3_2.pdf"))
    savefig(final_plot, joinpath(PLOTS_DIR, "figure3_combined.png"))
    savefig(final_plot, joinpath(PLOTS_DIR, "figure3_combined.pdf"))
    println("Figure 3 saved as PNG and PDF in '$PLOTS_DIR'.")
end

"""
    generate_figure4()

Generates and saves the plot for Figure 4, which is a 2x3 grid of filled contour plots
showing the total control effect (ψ) for six different epidemiological scenarios.
For this figure, the effect of treatment during the latent stage is turned OFF.
"""
function generate_figure4()
    println("Generating Figure 4: Control Effect (ψ) Contour Grid...")

    # Define the ranges for the plot axes, as per the project outline.[1]
    p_range = 0.0:0.02:1.0       # Bed net coverage
    epsilon_range = 0.0:0.02:1.0 # Proportion of nets treated

    # --- Panel [1, 1]: Low Prevalence, Low Biting Rate ---
    plot_1_1 = begin
        params = setup_scenario_parameters(1.2, 1.0)
        params.latent_treatment_switch = 0.0 # Turn OFF latent treatment effect
        psi_grid = zeros(length(p_range), length(epsilon_range))

        for (ip, p_val) in enumerate(p_range)
            for (ie, epsilon_val) in enumerate(epsilon_range)
                params.p = p_val
                params.epsilon = epsilon_val
                psi_grid[ip, ie] = calculate_psi(params)
            end
        end

        contourf(
            epsilon_range, p_range, psi_grid,
            title="Low Prevalence, Low Biting Rate",
            xlabel=L"\epsilon" * " (Treatment Proportion)",
            ylabel=L"p" * " (Bed Net Coverage)",
            colorbar_title=L"\psi",
            clims=(0, 1), c=:viridis
        )
    end
    println("Panel (1,1) done... ")

    # --- Panel [2, 3]: Low Prevalence, High Biting Rate ---
    plot_1_2 = begin
        params = setup_scenario_parameters(1.2, 2.0)
        params.latent_treatment_switch = 0.0 # Turn OFF latent treatment effect
        psi_grid = zeros(length(p_range), length(epsilon_range))

        for (ip, p_val) in enumerate(p_range)
            for (ie, epsilon_val) in enumerate(epsilon_range)
                params.p = p_val
                params.epsilon = epsilon_val
                psi_grid[ip, ie] = calculate_psi(params)
            end
        end

        contourf(
            epsilon_range, p_range, psi_grid,
            title="Low Prevalence, High Biting Rate",
            xlabel=L"\epsilon" * " (Treatment Proportion)",
            ylabel=L"p" * " (Bed Net Coverage)",
            colorbar_title=L"\psi",
            clims=(0, 1), c=:viridis
        )
    end
    println("Panel (1,2) done... ")

    # --- Panel [3, 2]: Moderate Prevalence, Low Biting Rate ---
    plot_2_1 = begin
        params = setup_scenario_parameters(5.0, 1.0)
        params.latent_treatment_switch = 0.0 # Turn OFF latent treatment effect
        psi_grid = zeros(length(p_range), length(epsilon_range))

        for (ip, p_val) in enumerate(p_range)
            for (ie, epsilon_val) in enumerate(epsilon_range)
                params.p = p_val
                params.epsilon = epsilon_val
                psi_grid[ip, ie] = calculate_psi(params)
            end
        end

        contourf(
            epsilon_range, p_range, psi_grid,
            title="Moderate Prevalence, Low Biting Rate",
            xlabel=L"\epsilon" * " (Treatment Proportion)",
            ylabel=L"p" * " (Bed Net Coverage)",
            colorbar_title=L"\psi",
            clims=(0, 1), c=:viridis
        )
    end
    println("Panel (2,1) done... ")

    # --- Panel [2, 2]: Moderate Prevalence, High Biting Rate ---
    plot_2_2 = begin
        params = setup_scenario_parameters(5.0, 2.0)
        params.latent_treatment_switch = 0.0 # Turn OFF latent treatment effect
        psi_grid = zeros(length(p_range), length(epsilon_range))

        for (ip, p_val) in enumerate(p_range)
            for (ie, epsilon_val) in enumerate(epsilon_range)
                params.p = p_val
                params.epsilon = epsilon_val
                psi_grid[ip, ie] = calculate_psi(params)
            end
        end

        contourf(
            epsilon_range, p_range, psi_grid,
            title="Moderate Prevalence, High Biting Rate",
            xlabel=L"\epsilon" * " (Treatment Proportion)",
            ylabel=L"p" * " (Bed Net Coverage)",
            colorbar_title=L"\psi",
            clims=(0, 1), c=:viridis
        )
    end
    println("Panel (2,2) done... ")

    # --- Panel [4, 2]: High Prevalence, Low Biting Rate ---
    plot_3_1 = begin
        params = setup_scenario_parameters(10.0, 1.0)
        params.latent_treatment_switch = 0.0 # Turn OFF latent treatment effect
        psi_grid = zeros(length(p_range), length(epsilon_range))

        for (ip, p_val) in enumerate(p_range)
            for (ie, epsilon_val) in enumerate(epsilon_range)
                params.p = p_val
                params.epsilon = epsilon_val
                psi_grid[ip, ie] = calculate_psi(params)
            end
        end

        contourf(
            epsilon_range, p_range, psi_grid,
            title="High Prevalence, Low Biting Rate",
            xlabel=L"\epsilon" * " (Treatment Proportion)",
            ylabel=L"p" * " (Bed Net Coverage)",
            colorbar_title=L"\psi",
            clims=(0, 1), c=:viridis
        )
    end
    println("Panel (3,1) done... ")

    # --- Panel [4, 3]: High Prevalence, High Biting Rate ---
    plot_3_2 = begin
        params = setup_scenario_parameters(10.0, 2.0)
        params.latent_treatment_switch = 0.0 # Turn OFF latent treatment effect
        psi_grid = zeros(length(p_range), length(epsilon_range))

        for (ip, p_val) in enumerate(p_range)
            for (ie, epsilon_val) in enumerate(epsilon_range)
                params.p = p_val
                params.epsilon = epsilon_val
                psi_grid[ip, ie] = calculate_psi(params)
            end
        end

        contourf(
            epsilon_range, p_range, psi_grid,
            title="High Prevalence, High Biting Rate",
            xlabel=L"\epsilon" * " (Treatment Proportion)",
            ylabel=L"p" * " (Bed Net Coverage)",
            colorbar_title=L"\psi",
            clims=(0, 1), c=:viridis
        )
    end
    println("Panel (3,2) done... ")

    # --- Combine all subplots into a 3x2 grid layout ---
    final_plot = plot(
        plot_1_1, plot_1_2,
        plot_2_1, plot_2_2,
        plot_3_1, plot_3_2,
        layout=(3, 2),
        size=(1400, 1500)
    )

    # --- Save all plots ---
    savefig(plot_1_1, joinpath(PLOTS_DIR, "figure4_panel_1_1.png"))
    savefig(plot_1_2, joinpath(PLOTS_DIR, "figure4_panel_1_2.png"))
    savefig(plot_2_1, joinpath(PLOTS_DIR, "figure4_panel_2_1.png"))
    savefig(plot_2_2, joinpath(PLOTS_DIR, "figure4_panel_2_2.png"))
    savefig(plot_3_1, joinpath(PLOTS_DIR, "figure4_panel_3_1.png"))
    savefig(plot_3_2, joinpath(PLOTS_DIR, "figure4_panel_3_2.png"))
    savefig(final_plot, joinpath(PLOTS_DIR, "figure4_combined.png"))
    savefig(final_plot, joinpath(PLOTS_DIR, "figure4_combined.pdf"))
    println("Figure 4 saved as PNG and PDF in '$PLOTS_DIR'.")
end

# Call the function to generate Figure 3
generate_figure3()

# Call the function to generate Figure 4
generate_figure4()