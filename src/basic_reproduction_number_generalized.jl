module BasicReproductionNumberGeneralized

export Parameters, compute_dfe_generalized, compute_R0_generalized, compute_endemic_equilibrium_generalized

using LinearAlgebra

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
    N_LM::Int          # number of latent stages for untreated mosquitoes
    N_LT::Int          # number of latent stages for treated mosquitoes
    sigma_M::Float64   # progression rate for untreated latent stages (day⁻¹)
    sigma_T::Float64   # progression rate for treated latent stages (day⁻¹)
end

# Compute the disease-free equilibrium (DFE) state vector.
# The state vector is ordered as follows:
# [ S_H, I_H, S_M, E₁ᴹ, ..., Eₙᴹ, I_M, S_T, E₁ᵀ, ..., Eₙᵀ, I_T ]
# where n = N_LM for untreated and n = N_LT for treated mosquitoes.
function compute_dfe_generalized(p::Parameters)
    # For humans, S_H = 1 and I_H = 0.
    # For mosquitoes, use the same 2x2 system:
    #     [ -(t_rate+g)   h     ] [ S*_M ] = [ -g ]
    #     [  t_rate    -(h+g)  ] [ S*_T ]   [  0 ]
    A = [-(p.t_rate + p.g) p.h;
        p.t_rate -(p.h + p.g)]
    b_vec = [-p.g, 0.0]
    S_M, S_T = A \ b_vec

    # Determine the dimension of the state vector:
    # Humans: 2 compartments.
    # Untreated mosquitoes: S_M (1) + N_LM latent compartments + I_M (1).
    # Treated mosquitoes: S_T (1) + N_LT latent compartments + I_T (1).
    n = 2 + (1 + p.N_LM + 1) + (1 + p.N_LT + 1)
    u0 = zeros(n)

    # Humans:
    u0[1] = 1.0   # S_H
    u0[2] = 0.0   # I_H

    # Untreated mosquitoes:
    u0[3] = S_M   # S_M
    # Latent untreated compartments and I_M remain 0.

    # Treated mosquitoes:
    # S_T is located after untreated compartments: index = 3 + N_LM + 1.
    idx_ST = 3 + p.N_LM + 1
    u0[idx_ST] = S_T
    # Latent treated compartments and I_T remain 0.

    return u0
end

# Compute the basic reproduction number (R₀) using the matrix method,
# generalized to an arbitrary number of latent stages.
function compute_R0_generalized(p::Parameters; analytical_method::Bool=true)
    # Obtain the DFE state.
    dfe = compute_dfe_generalized(p)
    S_H_star = dfe[1]
    S_M_star = dfe[3]
    idx_ST = 3 + p.N_LM + 1
    S_T_star = dfe[idx_ST]

    # Define infected compartments:
    # [ I_H, E₁ᴹ, E₂ᴹ, ..., Eₙᴹ, I_M, E₁ᵀ, E₂ᵀ, ..., Eₙᵀ, I_T ]
    n_inf = 1 + (p.N_LM + 1) + (p.N_LT + 1)
    F = zeros(n_inf, n_inf)
    V = zeros(n_inf, n_inf)

    # Define indices:
    idx_IH = 1
    idx_E_M_start = 2
    idx_I_M = idx_E_M_start + p.N_LM  # last untreated compartment (infectious)
    idx_E_T_start = idx_I_M + 1
    idx_I_T = idx_E_T_start + p.N_LT  # last treated compartment (infectious)

    # --- Construct the F matrix (new infections) ---
    # New infections in humans from infectious mosquitoes:
    F[idx_IH, idx_I_M] = p.m * p.a * p.b * S_H_star  # from untreated infectious
    F[idx_IH, idx_I_T] = p.m * p.a * p.b * S_H_star  # from treated infectious
    # New infections in mosquitoes (first latent stage) from infectious humans:
    F[idx_E_M_start, idx_IH] = p.a * p.c * S_M_star    # for untreated branch
    F[idx_E_T_start, idx_IH] = p.a * p.c * S_T_star    # for treated branch

    # --- Construct the V matrix (transitions) ---
    # Human compartment:
    V[idx_IH, idx_IH] = p.r

    # Untreated mosquito branch:
    # First latent compartment:
    V[idx_E_M_start, idx_E_M_start] = p.t_rate + p.sigma_M + p.g
    # Progression through latent stages:
    for i in 1:(p.N_LM-1)
        row = idx_E_M_start + i
        col = idx_E_M_start + i - 1
        V[row, col] = -p.sigma_M
        V[row, row] = p.sigma_M + p.g
    end
    # Transition from last latent untreated stage to infectious:
    V[idx_I_M, idx_I_M] = p.g
    V[idx_I_M, idx_E_M_start+p.N_LM-1] = -p.sigma_M

    # Extra coupling: some untreated mosquitoes get treated during latency.
    # This corresponds to an extra outflow from the first untreated latent compartment
    # to the first treated latent compartment at rate t_rate.
    V[idx_E_T_start, idx_E_M_start] = -p.t_rate

    # Treated mosquito branch:
    # For the first latent stage in treated branch:
    V[idx_E_T_start, idx_E_T_start] += p.sigma_T + p.g  # note the += because we already added the -t_rate above
    # Progression through treated latent stages:
    for i in 1:(p.N_LT-1)
        row = idx_E_T_start + i
        col = idx_E_T_start + i - 1
        V[row, col] = -p.sigma_T
        V[row, row] = p.sigma_T + p.g
    end
    # Transition from last treated latent stage to infectious treated compartment:
    V[idx_I_T, idx_I_T] = p.g
    V[idx_I_T, idx_E_T_start+p.N_LT-1] = -p.sigma_T

    # --- Next Generation Matrix Calculation ---
    V_inv = inv(V)
    NGM = F * V_inv
    eigenvalues = eigvals(NGM)
    R0 = maximum(real(eigenvalues))

    return R0
end

function compute_endemic_equilibrium_generalized(p::Parameters; tol=1e-8, max_iter=1000)
    # Define helper functions that compute the mosquito equilibrium components given I_H
    # Let D = a*c*I_H.
    function mosquito_components(I_H)
        D = p.a * p.c * I_H
        # Solve for S_M using the balance:
        # S_M = g / (D + t_rate + g - (h*t_rate/(D+h+g)) )
        S_M = p.g / (D + p.t_rate + p.g - (p.h * p.t_rate) / (D + p.h + p.g))
        S_T = p.t_rate * S_M / (D + p.h + p.g)
        E1M = (p.a * p.c * I_H * S_M) / (p.t_rate + p.sigma_M + p.g)
        # Propagate untreated latent stages:
        multiplier_M = p.sigma_M / (p.sigma_M + p.g)
        I_M = (p.sigma_M / p.g) * (multiplier_M)^(p.N_LM - 1) * E1M
        E1T = (p.a * p.c * I_H * S_T + p.t_rate * E1M) / (p.sigma_T + p.g)
        multiplier_T = p.sigma_T / (p.sigma_T + p.g)
        I_T = (p.sigma_T / p.g) * (multiplier_T)^(p.N_LT - 1) * E1T
        return S_M, S_T, E1M, I_M, E1T, I_T
    end

    # The consistency function that should be zero at equilibrium:
    # f(I_H) = I_H/(1-I_H) - m*a*b*(I_M+I_T)/r
    function consistency(I_H)
        _, _, _, I_M, _, I_T = mosquito_components(I_H)
        return I_H / (1 - I_H) - p.m * p.a * p.b * (I_M + I_T) / p.r
    end

    # Solve for I_H using the bisection method. We assume I_H is in (0,1).
    # We need to find a root of f(I_H) = 0.
    a_val, b_val = 1e-6, 0.999999
    f_a = consistency(a_val)
    f_b = consistency(b_val)
    if f_a * f_b > 0
        error("No sign change in consistency equation. Check parameter values.")
    end

    I_H_sol = 0.0
    for iter in 1:max_iter
        mid = (a_val + b_val) / 2
        f_mid = consistency(mid)
        if abs(f_mid) < tol
            I_H_sol = mid
            break
        end
        if f_mid * f_a < 0
            b_val = mid
            f_b = f_mid
        else
            a_val = mid
            f_a = f_mid
        end
        if iter == max_iter
            error("Bisection did not converge in compute_endemic_equilibrium_generalized")
        end
    end

    I_H = I_H_sol
    S_H = 1 - I_H
    S_M, S_T, E1M, I_M, E1T, I_T = mosquito_components(I_H)

    # Now compute the full latent stage compartments by propagating the ratios
    # for untreated mosquitoes:
    E_M = zeros(p.N_LM)
    E_M[1] = E1M
    for i in 2:p.N_LM
        E_M[i] = (p.sigma_M / (p.sigma_M + p.g)) * E_M[i-1]
    end
    # for treated mosquitoes:
    E_T = zeros(p.N_LT)
    E_T[1] = E1T
    for i in 2:p.N_LT
        E_T[i] = (p.sigma_T / (p.sigma_T + p.g)) * E_T[i-1]
    end

    # Assemble the state vector in the same order as the dfe vector
    # Order: [S_H, I_H, S_M, E_M (all stages), I_M, S_T, E_T (all stages), I_T]
    n = 2 + (1 + p.N_LM + 1) + (1 + p.N_LT + 1)
    eq = zeros(n)
    eq[1] = S_H
    eq[2] = I_H
    eq[3] = S_M
    for i in 1:p.N_LM
        eq[3+i] = E_M[i]
    end
    idx_IM = 3 + p.N_LM + 1
    eq[idx_IM] = I_M
    idx_ST = idx_IM + 1
    eq[idx_ST] = S_T
    for i in 1:p.N_LT
        eq[idx_ST+i] = E_T[i]
    end
    idx_IT = idx_ST + p.N_LT + 1
    eq[idx_IT] = I_T

    # @printf("Converged I_H = %.6f\n", I_H)
    return eq
end

end # module
