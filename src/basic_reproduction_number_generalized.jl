module BasicReproductionNumberGeneralized

export Parameters, compute_dfe_generalized, compute_R0_generalized, compute_endemic_equilibrium_generalized, generalized_ode!, compute_ee_numerical_generalized

using LinearAlgebra
using Roots
using Statistics

# Extended parameters: new fields for the number of latent stages and progression rates.
mutable struct Parameters
    a::Float64         # biting rate (day⁻¹)
    b::Float64         # human-to-mosquito transmission probability
    c::Float64         # mosquito-to-human transmission probability
    m::Float64         # mosquito-to-human ratio
    r::Float64         # human recovery rate (day⁻¹)
    g::Float64         # mosquito death rate (day⁻¹)
    h::Float64         # treatment waning rate (day⁻¹)
    t_rate::Float64    # treatment encounter rate (day⁻¹)
    L_NM::Int          # number of latent stages for untreated mosquitoes
    L_NT::Int          # number of latent stages for treated mosquitoes
    s_M::Float64       # total progression rate for untreated latent period (day⁻¹)
    s_T::Float64       # total progression rate for treated latent period (day⁻¹)
end

# Compute the disease-free equilibrium (DFE) state vector.
# The state vector is ordered as follows:
# [ S_H, I_H, S_M, E₁ᴹ, ..., E_L_NMᴹ, I_M, S_T, E₁ᵀ, ..., E_L_NTᵀ, I_T ]
function compute_dfe_generalized(p::Parameters)
    # For humans, S_H = 1 and I_H = 0.
    # For mosquitoes, solve the same 2x2 system:
    #     [ -(t_rate+g)   h     ] [ S*_M ] = [ -g ]
    #     [  t_rate    -(h+g)  ] [ S*_T ]   [  0 ]
    A = [-(p.t_rate + p.g) p.h;
        p.t_rate -(p.h + p.g)]
    b_vec = [-p.g, 0.0]
    S_M, S_T = A \ b_vec

    # Determine the dimension of the state vector:
    # Humans: 2 compartments.
    # Untreated mosquitoes: S_M (1) + L_NM latent compartments + I_M (1).
    # Treated mosquitoes: S_T (1) + L_NT latent compartments + I_T (1).
    n = 2 + (1 + p.L_NM + 1) + (1 + p.L_NT + 1)
    u0 = zeros(n)

    # Humans:
    u0[1] = 1.0   # S_H
    u0[2] = 0.0   # I_H

    # Untreated mosquitoes:
    u0[3] = S_M   # S_M; latent compartments and I_M remain 0.

    # Treated mosquitoes:
    # S_T is located after untreated compartments: index = 3 + L_NM + 1.
    idx_ST = 3 + p.L_NM + 1
    u0[idx_ST] = S_T  # latent compartments and I_T remain 0.

    return u0
end

# Compute the basic reproduction number (R₀) using the matrix method,
# generalized to an arbitrary number of latent stages.
function compute_R0_generalized(p::Parameters; analytical_method::Bool=true)
    # Obtain the DFE state.
    dfe = compute_dfe_generalized(p)
    S_H_star = dfe[1]
    S_M_star = dfe[3]
    idx_ST = 3 + p.L_NM + 1
    S_T_star = dfe[idx_ST]

    # Define infected compartments:
    # [ I_H, E₁ᴹ, E₂ᴹ, ..., E_L_NMᴹ, I_M, E₁ᵀ, E₂ᵀ, ..., E_L_NTᵀ, I_T ]
    n_inf = 1 + (p.L_NM + 1) + (p.L_NT + 1)
    F = zeros(n_inf, n_inf)
    V = zeros(n_inf, n_inf)

    # Define indices:
    idx_IH = 1
    idx_E_M_start = 2
    idx_I_M = idx_E_M_start + p.L_NM  # infectious untreated mosquito compartment
    idx_E_T_start = idx_I_M + 1
    idx_I_T = idx_E_T_start + p.L_NT  # infectious treated mosquito compartment

    # Per-stage progression rates:
    s_M_stage = p.s_M / p.L_NM
    s_T_stage = p.s_T / p.L_NT

    # --- Construct the F matrix (new infections) ---
    # New infections in humans from infectious mosquitoes:
    F[idx_IH, idx_I_M] = p.m * p.a * p.b * S_H_star  # from untreated infectious
    F[idx_IH, idx_I_T] = p.m * p.a * p.b * S_H_star  # from treated infectious
    # New infections in mosquitoes (first latent stage) from infectious humans:
    F[idx_E_M_start, idx_IH] = p.a * p.c * S_M_star    # untreated branch
    F[idx_E_T_start, idx_IH] = p.a * p.c * S_T_star    # treated branch

    # --- Construct the V matrix (transitions) ---
    # Human compartment:
    V[idx_IH, idx_IH] = p.r

    # Untreated mosquito branch:
    # First latent compartment:
    V[idx_E_M_start, idx_E_M_start] = p.t_rate + s_M_stage + p.g
    # Progression through untreated latent stages:
    for i in 1:(p.L_NM-1)
        row = idx_E_M_start + i
        col = idx_E_M_start + i - 1
        V[row, col] = -s_M_stage
        V[row, row] = s_M_stage + p.g
    end
    # Transition from last untreated latent stage to infectious:
    V[idx_I_M, idx_I_M] = p.g
    V[idx_I_M, idx_E_M_start+p.L_NM-1] = -s_M_stage

    # Extra coupling: untreated mosquitoes transferring to treated branch at rate t_rate.
    V[idx_E_T_start, idx_E_M_start] = -p.t_rate

    # Treated mosquito branch:
    # For the first latent stage:
    V[idx_E_T_start, idx_E_T_start] += s_T_stage + p.g   # note the += due to the -t_rate already present
    # Progression through treated latent stages:
    for i in 1:(p.L_NT-1)
        row = idx_E_T_start + i
        col = idx_E_T_start + i - 1
        V[row, col] = -s_T_stage
        V[row, row] = s_T_stage + p.g
    end
    # Transition from last treated latent stage to infectious:
    V[idx_I_T, idx_I_T] = p.g
    V[idx_I_T, idx_E_T_start+p.L_NT-1] = -s_T_stage

    # --- Next Generation Matrix Calculation ---
    V_inv = inv(V)
    NGM = F * V_inv
    eigenvalues = eigvals(NGM)
    R0 = maximum(real(eigenvalues))

    return R0
end

function compute_endemic_equilibrium_generalized(p::Parameters; tol=1e-8)
    # Helper function: compute mosquito equilibrium components given I_H.
    # Returns a tuple: (S_M, S_T, E₁ᴹ, I_M, E₁ᵀ, I_T)
    function mosquito_components(I_H)
        D = p.a * p.c * I_H
        S_M = p.g / (D + p.t_rate + p.g - (p.h * p.t_rate) / (D + p.h + p.g))
        S_T = p.t_rate * S_M / (D + p.h + p.g)
        s_M_stage = p.s_M / p.L_NM
        E1M = (p.a * p.c * I_H * S_M) / (p.t_rate + s_M_stage + p.g)
        multiplier_M = s_M_stage / (s_M_stage + p.g)
        I_M = (s_M_stage / p.g) * (multiplier_M)^(p.L_NM - 1) * E1M
        s_T_stage = p.s_T / p.L_NT
        E1T = (p.a * p.c * I_H * S_T + p.t_rate * E1M) / (s_T_stage + p.g)
        multiplier_T = s_T_stage / (s_T_stage + p.g)
        I_T = (s_T_stage / p.g) * (multiplier_T)^(p.L_NT - 1) * E1T
        return (S_M, S_T, E1M, I_M, E1T, I_T)
    end

    # Consistency function for the human equilibrium:
    # f(I_H) = I_H/(1-I_H) - m*a*b*(I_M+I_T)/r.
    function consistency(I_H)
        comps = mosquito_components(I_H)
        # Explicitly index the tuple:
        I_M = comps[4]
        I_T = comps[6]
        return I_H / (1 - I_H) - p.m * p.a * p.b * (I_M + I_T) / p.r
    end

    # Define the search interval:
    a, b = 1e-6, 1.0 - 1e-6
    fa = consistency(a)
    fb = consistency(b)
    if fa * fb < 0
        I_H_sol = find_zero(consistency, (a, b); atol=tol)
    else
        # Fall back to using an initial guess if the bracket does not work.
        I_H_sol = find_zero(consistency, 0.01; atol=tol)
    end

    I_H = I_H_sol
    S_H = 1 - I_H
    comps = mosquito_components(I_H)
    S_M, S_T, E1M, I_M, E1T, I_T = comps

    # Compute latent compartments for untreated mosquitoes.
    s_M_stage = p.s_M / p.L_NM
    E_M = zeros(p.L_NM)
    E_M[1] = E1M
    for i in 2:p.L_NM
        E_M[i] = (s_M_stage / (s_M_stage + p.g)) * E_M[i-1]
    end

    # Compute latent compartments for treated mosquitoes.
    s_T_stage = p.s_T / p.L_NT
    E_T = zeros(p.L_NT)
    E_T[1] = E1T
    for i in 2:p.L_NT
        E_T[i] = (s_T_stage / (s_T_stage + p.g)) * E_T[i-1]
    end

    # Assemble the full state vector.
    # Order: [S_H, I_H, S_M, E_M(1)...E_M(p.L_NM), I_M, S_T, E_T(1)...E_T(p.L_NT), I_T]
    n = 2 + (1 + p.L_NM + 1) + (1 + p.L_NT + 1)
    eq = zeros(n)
    eq[1] = S_H
    eq[2] = I_H
    eq[3] = S_M
    for i in 1:p.L_NM
        eq[3+i] = E_M[i]
    end
    idx_IM = 3 + p.L_NM + 1
    eq[idx_IM] = I_M
    idx_ST = idx_IM + 1
    eq[idx_ST] = S_T
    for i in 1:p.L_NT
        eq[idx_ST+i] = E_T[i]
    end
    idx_IT = idx_ST + p.L_NT + 1
    eq[idx_IT] = I_T

    return eq
end

using DifferentialEquations

# The generalized ODE system for the malaria model.
# State vector order:
# [ S_H, I_H, S_M, E₁ᴹ, …, Eₗᴹ, I_M, S_T, E₁ᵀ, …, Eₗᵀ, I_T ],
# where L = L_NM for untreated and L = L_NT for treated mosquitoes.
function generalized_ode!(du, u, p::Parameters, t)
    # Humans
    S_H = u[1]
    I_H = u[2]

    # Define indices:
    idx_SM = 3                                  # S_M index.
    idx_EM_start = idx_SM + 1                   # First untreated latent compartment.
    idx_IM = idx_EM_start + p.L_NM              # Infectious untreated (I_M).
    idx_ST = idx_IM + 1                        # S_T index.
    idx_ET_start = idx_ST + 1                   # First treated latent compartment.
    idx_IT = idx_ET_start + p.L_NT              # Infectious treated (I_T).

    # Retrieve mosquito compartments.
    S_M = u[idx_SM]
    S_T = u[idx_ST]
    I_M = u[idx_IM]
    I_T = u[idx_IT]

    # Per-stage progression rates.
    s_M_stage = p.s_M / p.L_NM
    s_T_stage = p.s_T / p.L_NT

    # -----------------------------
    # Human (SIS) dynamics:
    du[1] = -p.m * p.a * p.b * (I_M + I_T) * S_H + p.r * I_H
    du[2] = p.m * p.a * p.b * (I_M + I_T) * S_H - p.r * I_H

    # -----------------------------
    # Untreated mosquito dynamics:
    # Susceptible mosquitoes:
    du[idx_SM] = p.g + p.h * S_T - p.a * p.c * I_H * S_M - p.t_rate * S_M - p.g * S_M
    # Latent compartments:
    # E₁ᴹ:
    du[idx_EM_start] = p.a * p.c * I_H * S_M - (p.t_rate + s_M_stage + p.g) * u[idx_EM_start]
    # E₂ᴹ to E_Lᴹ:
    for i in 2:p.L_NM
        row = idx_EM_start + i - 1
        du[row] = s_M_stage * u[row-1] - (s_M_stage + p.g) * u[row]
    end
    # Infectious untreated mosquitoes:
    du[idx_IM] = s_M_stage * u[idx_EM_start+p.L_NM-1] - p.g * I_M

    # -----------------------------
    # Treated mosquito dynamics:
    # Susceptible mosquitoes:
    du[idx_ST] = p.t_rate * S_M - p.a * p.c * I_H * S_T - p.h * S_T - p.g * S_T
    # Latent compartments:
    # E₁ᵀ (first latent treated):
    du[idx_ET_start] = p.a * p.c * I_H * S_T + p.t_rate * u[idx_EM_start] - (s_T_stage + p.g) * u[idx_ET_start]
    # E₂ᵀ to E_Lᵀ:
    for i in 2:p.L_NT
        row = idx_ET_start + i - 1
        du[row] = s_T_stage * u[row-1] - (s_T_stage + p.g) * u[row]
    end
    # Infectious treated mosquitoes:
    du[idx_IT] = s_T_stage * u[idx_ET_start+p.L_NT-1] - p.g * I_T

    return du
end

# The ODE solver function that perturbs the DFE and integrates the system until convergence.
function compute_ee_numerical_generalized(p::Parameters; NLasts::Integer=100, δ::AbstractFloat=1e-5)
    # Initialize the state near the DFE (with a small seeding of infection)
    u0 = compute_dfe_generalized(p)
    u0[1] = 1 - 1e-3   # S_H slightly less than 1
    u0[2] = 1e-3       # I_H seeded with a small infection

    # Preallocate an array to store the last NLasts solutions.
    Lasts = Vector{Vector{Float64}}()

    # First, run a few iterations to initialize the Lasts array.
    for _ in 1:NLasts
        prob = ODEProblem(generalized_ode!, u0, (0.0, 1.0), p)
        sol = solve(prob, RK4(), dt=0.1, verbose=false)
        u0 = sol.u[end]
        push!(Lasts, deepcopy(u0))
    end

    # Iteratively integrate until the slope across the last NLasts is sufficiently small.
    while true
        # Combine the last solutions into a matrix (each column is a solution vector)
        M = hcat(Lasts...)  # size: (number of state vars) x NLasts

        # For each state variable (each row), compute the slope via simple linear regression.
        Slopes = map(eachrow(M)) do row
            N = length(row)
            X = collect(1:N)
            # Compute the slope: Cov(X, row) / Var(X)
            m_val = cov(X, row) / var(X)
            m_val
        end

        if all(abs.(Slopes) .< δ)
            break
        end

        # Integrate one more time from the last state.
        prob = ODEProblem(generalized_ode!, u0, (0.0, 1.0), p)
        sol = solve(prob, RK4(), dt=0.1, verbose=false)
        u0 = sol.u[end]
        popfirst!(Lasts)
        push!(Lasts, deepcopy(u0))
    end

    # Return the average of the last NLasts states as the endemic equilibrium.
    final_state = mean.(eachrow(hcat(Lasts...)))
    return final_state
end

end # module
