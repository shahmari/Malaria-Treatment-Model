 module MalariaTransmissionModel

using LinearAlgebra, NLsolve

export ModelParameters, malaria_ode!, calculate_R0, calculate_psi, calculate_S_M_star, calculate_pi_T, calculate_phi_T, calculate_R0_bar, compute_dfe, compute_R0_NGM, compute_endemic_equilibrium

"""
    ModelParameters

A struct to hold all fixed parameters of the malaria transmission model.
This approach ensures type stability and high performance for the ODE solver and other calculations.
"""
struct ModelParameters
    # Human Parameters
    r::Float64      # Human recovery rate (day^-1)

    # Mosquito Parameters
    a::Float64      # Biting rate (day^-1)
    m::Float64      # Mosquito-to-human ratio
    g::Float64      # Mosquito death rate (day^-1)
    L::Int          # Number of latent stages
    s_M::Float64    # Progression rate, untreated latency (day^-1)
    s_T::Float64    # Progression rate, treated latency (day^-1)

    # Transmission Parameters
    b::Float64      # Mosquito-to-human transmission probability
    c::Float64      # Human-to-mosquito transmission probability

    # Intervention Parameters
    h1::Float64     # Treatment waning rate 1 (day^-1)
    h2::Float64     # Treatment waning rate 2 (day^-1)
    theta::Float64  # Factor by which treatment increases incubation period
end

"""
    malaria_ode!(du, u, p_dynamic, t)

In-place function defining the 11-compartment ODE system for malaria transmission.
Designed for use with DifferentialEquations.jl solvers.
"""
function malaria_ode!(du, u, params::ModelParameters, t)
    # Unpack parameters for clarity
    # Note: p and epsilon are dynamic and passed in p_dynamic tuple
    p, epsilon = p_dynamic
    
    r, a, m, g, L, s_M, s_T, b, c, h1, h2 = 
        params.r, params.a, params.m, params.g, params.L, params.s_M, params.s_T, params.b, params.c, params.h1, params.h2

    # Unpack state variables for clarity
    S_H, I_H, S_M, E1_M, E2_M, I_M, S1_T, S2_T, E1_T, E2_T, I_T = u

    # Force of infection on humans
    foi_h = m * (1 - p) * a * b * (I_M + I_T)

    # Force of infection on mosquitoes
    foi_m = (1 - p) * a * c * I_H

    # Treatment acquisition rate for mosquitoes
    treat_rate = epsilon * p * a

    # System of Ordinary Differential Equations
    du[1] = -foi_h * S_H + r * I_H                                  # dS_H/dt
    du[2] = foi_h * S_H - r * I_H                                   # dI_H/dt
    du[3] = g + h2 * S2_T - foi_m * S_M - treat_rate * S_M - g * S_M # dS_M/dt
    du[4] = foi_m * S_M - (treat_rate + L * s_M + g) * E1_M         # dE1_M/dt
    du[5] = L * s_M * E1_M - (L * s_M + g) * E2_M                   # dE2_M/dt
    du[6] = L * s_M * E2_M - g * I_M                                # dI_M/dt
    du[7] = treat_rate * (S_M + S2_T) - (h1 + g) * S1_T             # dS1_T/dt
    du[8] = h1 * S1_T - (treat_rate + h2 + g) * S2_T                # dS2_T/dt
    du[9] = treat_rate * E1_M - (L * s_T + g) * E1_T                # dE1_T/dt
    du[10] = L * s_T * E1_T - (L * s_T + g) * E2_T                   # dE2_T/dt
    du[11] = L * s_T * E2_T - g * I_T                                # dI_T/dt
end

"""
    calculate_S_M_star(params, p, epsilon)

Calculates the proportion of susceptible, untreated mosquitoes at the 
disease-free equilibrium (DFE), S_M*.
Handles the case with and without waning treatment efficacy.
"""
function calculate_S_M_star(params::ModelParameters, p::Float64, epsilon::Float64)
    g, a, h1, h2 = params.g, params.a, params.h1, params.h2
    treat_rate = epsilon * p * a

    # General case with waning immunity (h1, h2 > 0)
    # Solved from the system of linear equations for DFE
    denominator = (h1 + g) * (treat_rate + h2 + g) - h1 * treat_rate
    if abs(denominator) < 1e-9 # Avoid division by zero
        # This case implies h1 and g are near zero, fallback to simpler model
        return g / (g + treat_rate)
    end
    
    num_S2_T = h1 * treat_rate * g
    den_S2_T = (treat_rate + g) * denominator - h1 * treat_rate * (g + h2)
    if abs(den_S2_T) < 1e-9
        S2_T_star = 0.0
    else
        S2_T_star = num_S2_T / den_S2_T
    end

    S_M_star = (g + h2 * S2_T_star) / (treat_rate + g)
    
    return S_M_star
end


"""
    calculate_pi_T(params, p, epsilon)

Calculates π_T, the probability that a mosquito acquires treatment during its latent period.
"""
function calculate_pi_T(params::ModelParameters, p::Float64, epsilon::Float64)
    a, L, s_M, g = params.a, params.L, params.s_M, params.g
    treat_rate = epsilon * p * a
    return treat_rate / (treat_rate + L * s_M + g)
end

"""
    calculate_phi_T(params)

Calculates ϕ_T, the factor by which treatment reduces the probability of a mosquito
surviving its (now extended) latent period.
"""
function calculate_phi_T(params::ModelParameters)
    L, s_M, s_T, g = params.L, params.s_M, params.s_T, params.g
    
    P_inv = ((L * s_M + g) / (L * s_M))^L
    surv_prob_treated = ((L * s_T) / (L * s_T + g))^L
    
    return P_inv * surv_prob_treated
end

"""
    calculate_psi(params, p, epsilon)

Calculates the total control effect, ψ = R₀ / R₀_bar.
This combines the three multiplicative factors related to the intervention.
"""
function calculate_psi(params::ModelParameters, p::Float64, epsilon::Float64)
    S_M_star = calculate_S_M_star(params, p, epsilon)
    pi_T = calculate_pi_T(params, p, epsilon)
    phi_T = calculate_phi_T(params)
    
    coverage_term = (1 - p)^2
    treatment_mix_term = (1 - pi_T) + phi_T * pi_T
    
    return treatment_mix_term * coverage_term * S_M_star
end

"""
    calculate_R0_bar(params)

Calculates the basic reproduction number without any control efforts, R₀_bar.
"""
function calculate_R0_bar(params)
    # Unpack parameters
    m, a, b, c, g, r, L, s_M = params.m, params.a, params.b, params.c, params.g, params.r, params.L, params.s_M

    # Probability of surviving the latent period (untreated)
    P = ((L * s_M) / (L * s_M + g))^L

    # Baseline R0 without any control efforts
    R0_bar = P * (m * a^2 * b * c) / (g * r)

    return R0_bar
end

"""
    calculate_R0(params, p, epsilon)

Calculates the basic reproduction number, R₀, using the symbolic Ross-Macdonald method.
"""
function calculate_R0(params::ModelParameters, p::Float64, epsilon::Float64)
    # Unpack parameters
    m, a, b, c, g, r, L, s_M = params.m, params.a, params.b, params.c, params.g, params.r, params.L, params.s_M
    
    # Probability of surviving the latent period (untreated)
    P = ((L * s_M) / (L * s_M + g))^L
    
    # Baseline R0 without any control efforts
    R0_bar = P * (m * a^2 * b * c) / (g * r)
    
    # Total control effect
    psi = calculate_psi(params, p, epsilon)
    
    # Final R0
    return psi * R0_bar
end

"""
    compute_dfe(params, p, epsilon)

Numerically computes the full 11-component state vector for the Disease-Free Equilibrium (DFE).
In the DFE, all infected compartments are zero, and the human population is fully susceptible.
"""
function compute_dfe(params::ModelParameters, p::Float64, epsilon::Float64)
    # At DFE, infected compartments are zero
    I_H, E1_M, E2_M, I_M, E1_T, E2_T, I_T = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    
    # All humans are susceptible
    S_H = 1.0

    # The susceptible mosquito compartments are at a steady state determined by treatment
    S_M_star = calculate_S_M_star(params, p, epsilon)
    
    # Unpack parameters needed to solve for the rest
    g, a, h1, h2 = params.g, params.a, params.h1, params.h2
    treat_rate = epsilon * p * a
    
    # Solve for S2_T and S1_T based on S_M_star
    # From S_M* = (g + h2*S2_T*)/(treat_rate + g) => S2_T* = (S_M* * (treat_rate + g) - g) / h2
    if abs(h2) > 1e-9
        S2_T_star = (S_M_star * (treat_rate + g) - g) / h2
    else
        S2_T_star = 0.0
    end

    # From S1_T* = treat_rate * (S_M* + S2_T*)/(h1 + g)
    if abs(h1 + g) > 1e-9
        S1_T_star = treat_rate * (S_M_star + S2_T_star) / (h1 + g)
    else
        S1_T_star = 0.0
    end

    # Return the full 11-element DFE state vector
    return
end


"""
    compute_R0_NGM(params, p, epsilon)

Computes the basic reproduction number, R₀, using the Next Generation Matrix (NGM) method.
This provides a numerical alternative to the symbolic `calculate_R0` function.
R₀ is the spectral radius (largest eigenvalue) of the matrix F * V⁻¹.
"""
function compute_R0_NGM(params::ModelParameters, p::Float64, epsilon::Float64)
    # Unpack parameters
    r, a, m, g, L, s_M, s_T, b, c = 
        params.r, params.a, params.m, params.g, params.L, params.s_M, params.s_T, params.b, params.c

    # R0 is calculated at the DFE
    S_M_star = calculate_S_M_star(params, p, epsilon)
    treat_rate = epsilon * p * a

    # Infected compartments are ordered as:
    
    # F matrix (new infections)
    F = zeros(7, 7)
    F[1, 2] = m * (1 - p) * a * b
    F[1, 3] = m * (1 - p) * a * b
    F[4, 1] = (1 - p) * a * c * S_M_star

    # V matrix (transitions between infected compartments)
    V = zeros(7, 7)
    V[1, 1] = r
    V[2, 2] = treat_rate + L * s_M + g
    V[5, 4] = -L * s_M
    V[3, 3] = L * s_M + g
    V[2, 5] = -L * s_M
    V[4, 4] = g
    V[6, 4] = -treat_rate
    V[5, 5] = L * s_T + g
    V[7, 6] = -L * s_T
    V[6, 6] = L * s_T + g
    V[3, 7] = -L * s_T
    V[7, 7] = g

    # Next Generation Matrix
    NGM = F * inv(V)

    # R0 is the spectral radius (magnitude of the largest eigenvalue)
    return maximum(abs.(eigvals(NGM)))
end


"""
    compute_endemic_equilibrium(params, p, epsilon; u0_guess=nothing)

Numerically solves for the endemic equilibrium state of the ODE system.
This is a non-trivial steady state where the disease persists in the population.
Uses a non-linear solver to find the roots of the ODE system.
"""
function compute_endemic_equilibrium(params::ModelParameters, p::Float64, epsilon::Float64; u0_guess=nothing)
    # Define the objective function for the root solver.
    # We want to find `u` such that `du/dt = 0`.
    function objective!(du, u)
        # Ensure state variables are non-negative, which is a physical constraint
        # and helps the solver.
        u[u.< 0].= 0
        
        # The total population must sum to 1. We enforce this by setting one
        # variable based on the others. Let's choose S_H.
        u[1] = 1.0 - sum(u[2:end])
        
        malaria_ode!(du, u, (p, epsilon), params)
    end

    # Provide a default initial guess if none is given.
    # A good guess is the DFE slightly perturbed towards an endemic state.
    if u0_guess === nothing
        u0 = compute_dfe(params, p, epsilon)
        u0[4] = 0.01 # Small number of infected humans
        u0[1] -= 0.01
    else
        u0 = u0_guess
    end

    # Use NLsolve to find the root of the objective function
    solution = nlsolve(objective!, u0)

    # Check for convergence
    if!converged(solution)
        @warn "The endemic equilibrium solver did not converge."
    end

    # Return the state vector at the endemic equilibrium
    return solution.zero
end

end # end of module MalariaTransmissionModel