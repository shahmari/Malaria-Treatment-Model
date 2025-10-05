module EffectiveReproductionNumber

export Parameters, compute_ee_numerical, malaria_ode!, compute_Re

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

function compute_ee_numerical(p::Parameters; NLasts::Integer=50, δ::AbstractFloat=1e-6)
    u0 = [0.8, 0.2, 0.3, 0.1, 0.1, 0.1, 0.2, 0.05, 0.05, 0.1]
    Lasts = Vector[]
    for _ ∈ 1:NLasts
        prob = ODEProblem(malaria_ode!, u0, (0.0, 1.0), p)
        sol = solve(prob, Rodas5(), abstol=1e-8)
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
        prob = ODEProblem(malaria_ode!, u0, (0.0, 2.0), p)
        sol = solve(prob, Rodas5(), abstol=1e-8)
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

function compute_Re(p::Parameters)
    # Compute EE
    ee = compute_ee_numerical(p)
    S_H, I_H = ee[1], ee[2]
    S_M, S_T = ee[3], ee[7]

    # Define infected compartments: [I_H, E1M, E2M, I_M, E1T, E2T, I_T]
    # Matrix F (new infections)
    F = zeros(7, 7)
    F[1, 4] = p.m * p.a * p.b * S_H  # I_H from I_M (at EE, S_H < 1)
    F[1, 7] = p.m * p.a * p.b * S_H  # I_H from I_T
    F[2, 1] = p.a * p.c * S_M        # E1M from I_H (S_M at EE)
    F[5, 1] = p.a * p.c * S_T        # E1T from I_H (S_T at EE)

    # Matrix V (transitions remain the same as parameters are constant)
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

    # Compute NGM and spectral radius
    V_inv = inv(V)
    NGM = F * V_inv
    R_eff = maximum(real(eigvals(NGM)))

    return R_eff
end

end