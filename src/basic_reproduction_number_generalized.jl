module BasicReproductionNumberGeneralized

export Parameters, compute_dfe_generalized, compute_R0_generalized

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
    L_NM::Int          # number of latent stages for untreated mosquitoes
    L_NT::Int          # number of latent stages for treated mosquitoes
    s_M::Float64   # total progression rate for untreated latent period (day⁻¹)
    s_T::Float64   # total progression rate for treated latent period (day⁻¹)
end

# Compute the disease-free equilibrium (DFE) state vector.
# The state vector is ordered as follows:
# [ S_H, I_H, S_M, E₁ᴹ, ..., Eₙᴹ, I_M, S_T, E₁ᵀ, ..., Eₙᵀ, I_T ]
# where n = L_NM for untreated and n = L_NT for treated mosquitoes.
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
    # Untreated mosquitoes: S_M (1) + L_NM latent compartments + I_M (1).
    # Treated mosquitoes: S_T (1) + L_NT latent compartments + I_T (1).
    n = 2 + (1 + p.L_NM + 1) + (1 + p.L_NT + 1)
    u0 = zeros(n)

    # Humans:
    u0[1] = 1.0   # S_H
    u0[2] = 0.0   # I_H

    # Untreated mosquitoes:
    u0[3] = S_M   # S_M
    # Latent untreated compartments and I_M remain 0.

    # Treated mosquitoes:
    # S_T is located after untreated compartments: index = 3 + L_NM + 1.
    idx_ST = 3 + p.L_NM + 1
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
    idx_ST = 3 + p.L_NM + 1
    S_T_star = dfe[idx_ST]

    # Define infected compartments:
    # [ I_H, E₁ᴹ, E₂ᴹ, ..., Eₙᴹ, I_M, E₁ᵀ, E₂ᵀ, ..., Eₙᵀ, I_T ]
    n_inf = 1 + (p.L_NM + 1) + (p.L_NT + 1)
    F = zeros(n_inf, n_inf)
    V = zeros(n_inf, n_inf)

    # Define indices:
    idx_IH = 1
    idx_E_M_start = 2
    idx_I_M = idx_E_M_start + p.L_NM  # last untreated compartment (infectious)
    idx_E_T_start = idx_I_M + 1
    idx_I_T = idx_E_T_start + p.L_NT  # last treated compartment (infectious)

    # Pre-calculate per-stage progression rates:
    s_M_stage = p.s_M / p.L_NM
    s_T_stage = p.s_T / p.L_NT

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

    # Extra coupling: some untreated mosquitoes get treated during latency.
    # This corresponds to an extra outflow from the first untreated latent compartment
    # to the first treated latent compartment at rate t_rate.
    V[idx_E_T_start, idx_E_M_start] = -p.t_rate

    # Treated mosquito branch:
    # For the first latent stage in treated branch:
    V[idx_E_T_start, idx_E_T_start] += s_T_stage + p.g  # note the += because we already added the -t_rate above
    # Progression through treated latent stages:
    for i in 1:(p.L_NT-1)
        row = idx_E_T_start + i
        col = idx_E_T_start + i - 1
        V[row, col] = -s_T_stage
        V[row, row] = s_T_stage + p.g
    end
    # Transition from last treated latent stage to infectious treated compartment:
    V[idx_I_T, idx_I_T] = p.g
    V[idx_I_T, idx_E_T_start+p.L_NT-1] = -s_T_stage

    # --- Next Generation Matrix Calculation ---
    V_inv = inv(V)
    NGM = F * V_inv
    eigenvalues = eigvals(NGM)
    R0 = maximum(real(eigenvalues))

    return R0
end

end # module
