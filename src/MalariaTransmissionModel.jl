module MalariaTransmissionModel

using LinearAlgebra, NLsolve

export ModelParameters, malaria_ode!, calculate_R0, calculate_psi, 
    calculate_pi_T, calculate_phi_T, calculate_R0_bar, compute_dfe, 
    compute_R0_NGM, compute_endemic_equilibrium, compute_dfe_numeric

"""
    ModelParameters

A struct to hold all parameters of the malaria transmission model.
"""
mutable struct ModelParameters
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

    # Control Parameters
    p::Float64      # Coverage
    epsilon::Float64 # Efficacy

    # NEW: Switch for latent-stage treatment
    latent_treatment_switch::Float64 # 1.0 for on (default), 0.0 for off
end

"""
    malaria_ode!(du, u, params, t)

In-place function defining the 11-compartment ODE system for malaria transmission.
"""
function malaria_ode!(du, u, params, t)
    # Unpack all parameters from the struct
    r, a, m, g, L, s_M, s_T, b, c, h1, h2, p, epsilon, latent_switch = 
     params.r, params.a, params.m, params.g, params.L, params.s_M, params.s_T, 
     params.b, params.c, params.h1, params.h2, params.p, params.epsilon, 
     params.latent_treatment_switch

    # Unpack state variables
    S_H, I_H, S_M, E1_M, E2_M, I_M, S1_T, S2_T, E1_T, E2_T, I_T = u

    # Rates
    foi_h = m * (1 - p) * a * b * (I_M + I_T)
    foi_m = (1 - p) * a * c * I_H
    treat_rate = epsilon * p * a

    # System of Ordinary Differential Equations
    du[1] = -foi_h * S_H + r * I_H                                  # dS_H/dt
    du[2] = foi_h * S_H - r * I_H                                   # dI_H/dt
    du[3] = g + h2 * S2_T - foi_m * S_M - treat_rate * S_M - g * S_M # dS_M/dt
    
    # MODIFIED: The outflow from E1_M due to treatment is now controlled by the switch
    du[4] = foi_m * S_M - (latent_switch * treat_rate + L * s_M + g) * E1_M # dE1_M/dt
    
    du[5] = L * s_M * E1_M - (L * s_M + g) * E2_M                   # dE2_M/dt
    du[6] = L * s_M * E2_M - g * I_M                                # dI_M/dt
    du[7] = treat_rate * (S_M + S2_T) - (h1 + g) * S1_T             # dS1_T/dt
    du[8] = h1 * S1_T - (treat_rate + h2 + g) * S2_T                # dS2_T/dt

    # MODIFIED: The inflow to E1_T from E1_M is now controlled by the switch
    du[9] = latent_switch * treat_rate * E1_M - (L * s_T + g) * E1_T # dE1_T/dt

    du[10] = L * s_T * E1_T - (L * s_T + g) * E2_T                   # dE2_T/dt
    du[11] = L * s_T * E2_T - g * I_T                                # dI_T/dt
end

"""
    calculate_pi_T(params)

Calculates π_T, the probability that a mosquito acquires treatment during its latent period.
"""
function calculate_pi_T(params::ModelParameters)
    a, L, s_M, g, p, epsilon, latent_switch = 
        params.a, params.L, params.s_M, params.g, params.p, params.epsilon, params.latent_treatment_switch
    
    treat_rate = epsilon * p * a
    
    # MODIFIED: The rate of acquiring treatment for a latent mosquito is controlled by the switch
    switched_treat_rate = latent_switch * treat_rate
    
    return switched_treat_rate / (switched_treat_rate + L * s_M + g)
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
    calculate_psi(params)

Calculates the total control effect, ψ = R₀ / R₀_bar.
"""
function calculate_psi(params::ModelParameters)
    S_M_star = compute_dfe(params).S_M
    pi_T = calculate_pi_T(params)
    phi_T = calculate_phi_T(params)
    p = params.p

    coverage_term = (1 - p)^2
    treatment_mix_term = (1 - pi_T) + phi_T * pi_T
    
    return treatment_mix_term * coverage_term * S_M_star
end

"""
    calculate_R0_bar(params)

Calculates the basic reproduction number without any control efforts, R₀_bar.
"""
function calculate_R0_bar(params::ModelParameters)
    m, a, b, c, g, r, L, s_M = params.m, params.a, params.b, params.c, params.g, params.r, params.L, params.s_M
    P = ((L * s_M) / (L * s_M + g))^L
    return P * (m * a^2 * b * c) / (g * r)
end

"""
    calculate_R0(params)

Calculates the basic reproduction number, R₀, using the symbolic Ross-Macdonald method.
"""
function calculate_R0(params::ModelParameters)
    R0_bar = calculate_R0_bar(params)
    psi = calculate_psi(params)
    return psi * R0_bar
end

"""
    compute_dfe(params)

Computes the full 11-component state vector for the Disease-Free Equilibrium (DFE).
"""
function compute_dfe(params::ModelParameters)
    a, g, h1, h2, p, epsilon = params.a, params.g, params.h1, params.h2, params.p, params.epsilon
    treat_rate = epsilon * p * a

    S1_T_star = treat_rate / (treat_rate + h1 + g)
    S2_T_star = (h1 / (treat_rate + h2 + g)) * S1_T_star
    S_M_star_val = 1.0 - S1_T_star - S2_T_star

    return (S_H=1.0, I_H=0.0, S_M=S_M_star_val, E1_M=0.0, E2_M=0.0, I_M=0.0, 
            S1_T=S1_T_star, S2_T=S2_T_star, E1_T=0.0, E2_T=0.0, I_T=0.0)
end

"""
    compute_dfe_numeric(params)

Numerically computes the disease-free equilibrium (DFE) state.
"""
function compute_dfe_numeric(params::ModelParameters)
    S_H = 1.0
    I_H = 0.0

    function f!(F, vars)
        S_M, S2_T = vars
        a, g, h1, h2, p, epsilon = params.a, params.g, params.h1, params.h2, params.p, params.epsilon
        treat_rate = epsilon * p * a
        S1_T = treat_rate * (S_M + S2_T) / (h1 + g)
        F[1] = g + h2 * S2_T - ((treat_rate + g) * S_M)
        F[2] = h1 * S1_T - ((treat_rate * S2_T + h2 + g) * S2_T)
    end

    guess = [params.g / (params.g + params.epsilon * params.p * params.a), 0.0]
    sol = nlsolve(f!, guess)
    S_M, S2_T = sol.zero
    S1_T = params.epsilon * params.p * params.a * (S_M + S2_T) / (params.h1 + params.g)

    return (S_H=S_H, I_H=I_H, S_M=S_M, E1_M=0.0, E2_M=0.0, I_M=0.0,
            S1_T=S1_T, S2_T=S2_T, E1_T=0.0, E2_T=0.0, I_T=0.0)
end

"""
    compute_R0_NGM(params)

Computes R₀ using the Next Generation Matrix (NGM) method.
"""
function compute_R0_NGM(params::ModelParameters)
    r, a, m, g, L, s_M, s_T, b, c, p, epsilon, latent_switch = 
        params.r, params.a, params.m, params.g, params.L, params.s_M, params.s_T, 
        params.b, params.c, params.p, params.epsilon, params.latent_treatment_switch
    
    S_M_star = compute_dfe(params).S_M
    treat_rate = epsilon * p * a

    F = zeros(7, 7)
    F[1, 4] = m * (1 - p) * a * b
    F[1, 7] = m * (1 - p) * a * b
    F[2, 1] = (1 - p) * a * c * S_M_star

    V = zeros(7, 7)
    V[1, 1] = r
    # MODIFIED: Outflow from E1_M and inflow to E1_T are controlled by the switch
    V[2, 2] = latent_switch * treat_rate + L * s_M + g
    V[3, 2] = -L * s_M
    V[3, 3] = L * s_M + g
    V[4, 3] = -L * s_M
    V[4, 4] = g
    V[5, 2] = -latent_switch * treat_rate
    V[5, 5] = L * s_T + g
    V[6, 5] = -L * s_T
    V[6, 6] = L * s_T + g
    V[7, 6] = -L * s_T
    V[7, 7] = g

    NGM = F * inv(V)
    return maximum(abs.(eigvals(NGM)))
end

"""
    compute_endemic_equilibrium(params; u0_guess=nothing)

Numerically solves for the endemic equilibrium state of the ODE system.
"""
function compute_endemic_equilibrium(params::ModelParameters; u0_guess=nothing)
    a, b, c, m = params.a, params.b, params.c, params.m
    r, g, h1, h2 = params.r, params.g, params.h1, params.h2
    L, s_M, s_T = params.L, params.s_M, params.s_T
    p, epsilon, latent_switch = params.p, params.epsilon, params.latent_treatment_switch

    λ1M = L * s_M
    λ2M = L * s_M
    λ1T = L * s_T
    λ2T = L * s_T
    τ = epsilon * p * a

    function resid!(F, x)
        SM, S2T, IH = x
        S1T = τ * (SM + S2T) / (h1 + g)

        # MODIFIED: Switched treatment rate applied to latent compartments
        E1M = (1 - p) * a * c * IH * SM / (latent_switch * τ + λ1M + g)
        E2M = λ1M * E1M / (λ2M + g)
        IM  = λ2M * E2M / g

        E1T = (latent_switch * τ * E1M) / (λ1T + g)
        E2T = λ1T * E1T / (λ2T + g)
        IT  = λ2T * E2T / g

        F[1] = g + h2 * S2T - ((1 - p) * a * c * IH + τ + g) * SM
        F[2] = h1 * S1T - ((τ * S2T + h2 + g) * S2T)
        F[3] = m * (1 - p) * a * b * (IM + IT) * (1 - IH) - r * IH
    end

    dfe = compute_dfe(params)
    guess = [dfe.S_M, 0.0, 1.0]

    sol = nlsolve(resid!, guess; ftol=1e-10, xtol=1e-10)
    if !converged(sol)
        @warn "Endemic-equilibrium solver failed to converge"
    end

    SM, S2_T, IH = sol.zero
    SH = 1.0 - IH

    S1T = τ * (SM + S2_T) / (h1 + g)
    E1M = (1 - p) * a * c * IH * SM / (latent_switch * τ + λ1M + g)
    E2M = λ1M * E1M / (λ2M + g)
    IM  = λ2M * E2M / g
    E1T = (latent_switch * τ * E1M) / (λ1T + g)
    E2T = λ1T * E1T / (λ2T + g)
    IT  = λ2T * E2T / g

    return (S_H=SH, I_H=IH, S_M=SM, E1_M=E1M, E2_M=E2M, I_M=IM,
            S1_T=S1T, S2_T=S2_T, E1_T=E1T, E2_T=E2T, I_T=IT)
end

end # end of module MalariaTransmissionModel