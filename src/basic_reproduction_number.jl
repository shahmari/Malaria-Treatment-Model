module BasicReproductionNumber

export Parameters, compute_dfe_matrix, compute_dfe_numerical, malaria_ode!, compute_R0, compute_endemic_equilibrium, compute_R0_symbolic

using LinearAlgebra, DifferentialEquations, Statistics

# Define a mutable struct for parameters
mutable struct Parameters
    a::Float64
    b::Float64
    c::Float64
    m::Float64
    r::Float64
    g::Float64
    h::Float64
    t_rate::Float64
    s1M::Float64
    s2M::Float64
    s1T::Float64
    s2T::Float64
end

function compute_dfe_matrix(p::Parameters)
    g, h, t_rate = p.g, p.h, p.t_rate

    # Solve for S_M and S_T at DFE (I_H = 0)
    A = [-(t_rate + g) h;
        t_rate -(h + g)]
    b_vec = [-g, 0.0]
    S_M, S_T = A \ b_vec

    return [1.0, 0.0, S_M, 0.0, 0.0, 0.0, S_T, 0.0, 0.0, 0.0]
end

function compute_dfe_numerical(p::Parameters; NLasts::Integer=100, δ::AbstractFloat=1e-5)
    u0 = [1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0]
    Lasts = Vector[]
    for _ ∈ 1:NLasts
        prob = ODEProblem(malaria_ode!, u0, (0.0, 1.0), p)
        sol = solve(prob, RK4(), verbose=false)
        u0 = sol.u[end]
        push!(Lasts, deepcopy(u0))
    end
    while true
        Slope = map(eachrow(hcat(Lasts...))) do Vec
            A = [hcat(1:NLasts) reshape(ones(NLasts), NLasts, 1)]
            b = reshape(Vec, NLasts, 1)
            x = A \ b
            x[1]
        end
        if all(abs.(Slope) .< δ)
            break
        end
        prob = ODEProblem(malaria_ode!, u0, (0.0, 1.0), p)
        sol = solve(prob, RK4(), verbose=false)
        u0 = sol.u[end]
        popfirst!(Lasts)
        push!(Lasts, deepcopy(u0))
    end
    return mean.(eachrow(hcat(Lasts...)))
end

function malaria_ode!(du, u, p::Parameters, t)
    # Unpack state variables
    S_H, I_H, S_M, E1M, E2M, I_M, S_T, E1T, E2T, I_T = u

    # Human dynamics
    du[1] = -p.m * p.a * p.b * (I_M + I_T) * S_H + p.r * I_H  # dS_H/dt
    du[2] = p.m * p.a * p.b * (I_M + I_T) * S_H - p.r * I_H   # dI_H/dt

    # Untreated mosquito dynamics
    du[3] = p.g + p.h * S_T - p.a * p.c * I_H * S_M - p.t_rate * S_M - p.g * S_M  # dS_M/dt
    du[4] = p.a * p.c * I_H * S_M - p.t_rate * E1M - p.s1M * E1M - p.g * E1M      # dE1M/dt
    du[5] = p.s1M * E1M - p.s2M * E2M - p.g * E2M                                 # dE2M/dt
    du[6] = p.s2M * E2M - p.g * I_M                                               # dI_M/dt

    # Treated mosquito dynamics
    du[7] = p.t_rate * S_M - p.a * p.c * I_H * S_T - p.h * S_T - p.g * S_T        # dS_T/dt
    du[8] = p.a * p.c * I_H * S_T + p.t_rate * E1M - p.s1T * E1T - p.g * E1T      # dE1T/dt
    du[9] = p.s1T * E1T - p.s2T * E2T - p.g * E2T                                 # dE2T/dt
    du[10] = p.s2T * E2T - p.g * I_T                                              # dI_T/dt

    return du
end

function compute_R0(p::Parameters; analytical_method::Bool=true)
    # Compute DFE
    dfe = analytical_method ? compute_dfe_matrix(p) : compute_dfe_numerical(p)
    S_H_star, I_H_star = dfe[1], dfe[2]
    S_M_star, S_T_star = dfe[3], dfe[7]

    # Ensure human population is normalized (S_H + I_H = 1)
    # @assert S_H_star ≈ 1.0 && I_H_star ≈ 0.0 "DFE violates S_H + I_H = 1"

    # Define infected compartments: [I_H, E1M, E2M, I_M, E1T, E2T, I_T]
    # Matrix F (new infections)
    F = zeros(7, 7)

    # Mosquito-to-human transmission (I_H infected by I_M and I_T)
    F[1, 4] = p.m * p.a * p.b * S_H_star  # From I_M (S_H = 1 at DFE)
    F[1, 7] = p.m * p.a * p.b * S_H_star  # From I_T

    # Human-to-mosquito transmission (E1M and E1T infected by I_H)
    F[2, 1] = p.a * p.c * S_M_star  # E1M from I_H
    F[5, 1] = p.a * p.c * S_T_star  # E1T from I_H

    # Matrix V (transitions)
    V = diagm([
        p.r,                   # I_H recovery
        p.t_rate + p.s1M + p.g, # E1M outflow (now includes t_rate)
        p.s2M + p.g,           # E2M outflow
        p.g,                    # I_M death
        p.s1T + p.g,            # E1T outflow
        p.s2T + p.g,            # E2T outflow
        p.g                     # I_T death
    ])

    # Off-diagonal transitions
    V[3, 2] = -p.s1M  # E1M -> E2M
    V[4, 3] = -p.s2M  # E2M -> I_M
    V[5, 2] = -p.t_rate  # E1M -> E1T (NEW CORRECTION)
    V[6, 5] = -p.s1T  # E1T -> E2T
    V[7, 6] = -p.s2T  # E2T -> I_T

    # Compute NGM and R0
    V_inv = inv(V)
    NGM = F * V_inv
    R0 = maximum(real(eigvals(NGM)))  # Use real part for stability

    return R0
end

function compute_endemic_equilibrium(p::Parameters)
    # Compute endemic equilibrium values
    SM = (p.g * (p.a * p.c + p.h + p.g)) / ((p.a * p.c + p.t_rate + p.g) * (p.a * p.c + p.h + p.g) - p.h * p.t_rate)
    ST = (p.t_rate / (p.a * p.c + p.h + p.g)) * SM

    IH = (p.m * p.a * p.b * (p.s2M / p.g) * (p.s1M / (p.s2M + p.g)) * (p.a * p.c * SM / (p.t_rate + p.s1M + p.g)) +
          p.m * p.a * p.b * (p.s2T / p.g) * (p.s1T / (p.s2T + p.g)) * (p.a * p.c * ST / (p.s1T + p.g))) /
         (p.r + p.m * p.a * p.b * (p.s2M / p.g) * (p.s1M / (p.s2M + p.g)) * (p.a * p.c * SM / (p.t_rate + p.s1M + p.g)) +
          p.m * p.a * p.b * (p.s2T / p.g) * (p.s1T / (p.s2T + p.g)) * (p.a * p.c * ST / (p.s1T + p.g)))
    SH = 1 - IH

    E1M = (p.a * p.c * IH * SM) / (p.t_rate + p.s1M + p.g)
    E2M = (p.s1M * E1M) / (p.s2M + p.g)
    IM = (p.s2M * E2M) / p.g

    E1T = (p.a * p.c * IH * ST + p.t_rate * E1M) / (p.s1T + p.g)
    E2T = (p.s1T * E1T) / (p.s2T + p.g)
    IT = (p.s2T * E2T) / p.g

    return [SH, IH, SM, E1M, E2M, IM, ST, E1T, E2T, IT]
end

function compute_R0_symbolic(p::Parameters)
    term1 = ((p.h + p.g) * p.s1M * p.s2M) / ((p.t_rate + p.h + p.g) * (p.t_rate + p.s1M + p.g) * (p.s2M + p.g))
    term2 = (p.t_rate * p.s1T * p.s2T) / ((p.t_rate + p.h + p.g) * (p.s1T + p.g) * (p.s2T + p.g))
    term3 = ((p.h + p.g) * p.t_rate * p.s1T * p.s2T) / ((p.t_rate + p.h + p.g) * (p.t_rate + p.s1M + p.g) * (p.s1T + p.g) * (p.s2T + p.g))

    R0 = (p.m * p.a^2 * p.b * p.c) / (p.r * p.g) * (term1 + term2 + term3)

    return R0
end


end