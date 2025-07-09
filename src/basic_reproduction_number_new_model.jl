module MalariaTreatmentModel

export Parameters,
    compute_dfe,
    compute_dfe_closed_form,
    compute_R0,
    compute_R0_closed_form,
    compute_endemic_equilibrium,
    compute_ee_numerical,
    malaria_ode!

using LinearAlgebra, NLsolve, DifferentialEquations, Statistics

"""
Parameters of the malaria model with mosquito treatment (2 latent stages).

Fields:
  a    :: Float64  # biting rate (day⁻¹)
  b    :: Float64  # human-to-mosquito transmission probability
  c    :: Float64  # mosquito-to-human transmission probability
  m    :: Float64  # mosquito-to-human ratio
  r    :: Float64  # human recovery rate (day⁻¹)
  g    :: Float64  # mosquito death rate (day⁻¹)
  h1   :: Float64  # treated susceptible transition rate (S1,T -> S2,T)
  h2   :: Float64  # treatment waning rate (S2,T -> SM)
  ε    :: Float64  # efficacy factor of treatment on mosquito biting
  p    :: Float64  # treatment coverage (fraction of bites affecting treatment)
  s1M  :: Float64  # progression rate untreated latent 1 -> 2 (day⁻¹)
  s2M  :: Float64  # progression rate untreated latent 2 -> I_M (day⁻¹)
  s1T  :: Float64  # progression rate treated latent 1 -> 2 (day⁻¹)
  s2T  :: Float64  # progression rate treated latent 2 -> I_T (day⁻¹)
"""
mutable struct Parameters
    a::Float64
    b::Float64
    c::Float64
    m::Float64
    r::Float64
    g::Float64
    h1::Float64
    h2::Float64
    ε::Float64
    p::Float64
    s1M::Float64
    s2M::Float64
    s1T::Float64
    s2T::Float64
end

"""
Compute the disease-free equilibrium (DFE) state numerically by solving the steady-state subsystem.
Returns a NamedTuple with fields:
  :SH, :IH,
  :SM, :E1M, :E2M, :IM,
  :S1T, :S2T, :E1T, :E2T, :IT
"""
function compute_dfe(par::Parameters)
    SH = 1.0
    IH = 0.0
    function f!(F, vars)
        SM, S2T = vars
        a = par.a; b = par.b; c = par.c; m = par.m;
        r = par.r; g = par.g; h1 = par.h1; h2 = par.h2;
        ε = par.ε; p = par.p; s1M = par.s1M;
        s2M = par.s2M; s1T = par.s1T; s2T = par.s2T
        S1T = ε * p * a * (SM + S2T) / (h1 + g)
        F[1] = g + h2 * S2T - ((ε * p * a + g) * SM)
        F[2] = h1 * S1T - ((ε * p * a * S2T + h2 + g) * S2T)
    end
    guess = [par.g / (par.g + par.ε * par.p * par.a), 0.0]
    sol = nlsolve(f!, guess)
    SM, S2T = sol.zero
    S1T = par.ε * par.p * par.a * (SM + S2T) / (par.h1 + par.g)
    return (SH=SH, IH=IH,
        SM=SM, E1M=0.0, E2M=0.0, IM=0.0,
        S1T=S1T, S2T=S2T,
        E1T=0.0, E2T=0.0, IT=0.0)
end

"""
Compute the disease-free equilibrium (DFE) state in closed form using algebraic expressions:
  S1T* = ε p a / (ε p a + h1 + g)
  S2T* = (h1 / (ε p a + h2 + g)) * S1T*
  SM*  = 1 - S1T* - S2T*
"""
function compute_dfe_closed_form(par::Parameters)
    SH = 1.0
    IH = 0.0
    a = par.a; b = par.b; c = par.c; m = par.m;
    r = par.r; g = par.g; h1 = par.h1; h2 = par.h2;
    ε = par.ε; p = par.p; s1M = par.s1M;
    s2M = par.s2M; s1T = par.s1T; s2T = par.s2T
    # Closed-form treated compartments
    S1T = ε * p * a / (ε * p * a + h1 + g)
    S2T = (h1 / (ε * p * a + h2 + g)) * S1T
    SM = 1.0 - S1T - S2T
    return (SH=SH, IH=IH,
        SM=SM, E1M=0.0, E2M=0.0, IM=0.0,
        S1T=S1T, S2T=S2T,
        E1T=0.0, E2T=0.0, IT=0.0)
end

"""
Compute the basic reproduction number R₀ via the next-generation method.
"""
function compute_R0(par::Parameters)
    dfe = compute_dfe_closed_form(par)
    SH, SM = dfe.SH, dfe.SM
    a = par.a; b = par.b; c = par.c; m = par.m;
    r = par.r; g = par.g; h1 = par.h1; h2 = par.h2;
    ε = par.ε; p = par.p; s1M = par.s1M;
    s2M = par.s2M; s1T = par.s1T; s2T = par.s2T
    F = zeros(7, 7)
    V = zeros(7, 7)
    F[1, 4] = m * (1 - p) * a * b * SH
    F[1, 7] = m * (1 - p) * a * b * SH
    F[2, 1] = (1 - p) * a * c * SM
    V[1, 1] = r
    V[2, 2] = ε * p * a + s1M + g
    V[3, 3] = s2M + g
    V[3, 2] = -s1M
    V[4, 4] = g
    V[4, 3] = -s2M
    V[5, 5] = s1T + g
    V[5, 2] = -ε * p * a
    V[6, 6] = s2T + g
    V[6, 5] = -s1T
    V[7, 7] = g
    V[7, 6] = -s2T
    K = F * (V \ I)
    return maximum(real(eigvals(K)))
end

function compute_R0_closed_form(par::Parameters)
    # Closed-form R0 from the report
    a = par.a; b = par.b; c = par.c; m = par.m;
    r = par.r; g = par.g; h1 = par.h1; h2 = par.h2;
    ε = par.ε; p = par.p; s1M = par.s1M;
    s2M = par.s2M; s1T = par.s1T; s2T = par.s2T
    # Survival probability through 2-stage untreated latency
    P = (s1M / (s1M + g)) * (s2M / (s2M + g))
    # Baseline R0 without treatment
    R0_bar = P * m * a^2 * b * c / (g * r)
    # Probability of acquiring treatment during latency (two stages): rate εpa over total exit rate
    pi_T = ε * p * a / (ε * p * a + s1M + s2M + g)
    # Survival probability through 2-stage treated latency
    P_T = (s1T / (s1T + g)) * (s2T / (s2T + g))
    # Factor phi_T = relative survival under treatment vs untreated
    phi_T = P_T / P
    # DFE susceptible untreated mosquito fraction
    SM = compute_dfe_closed_form(par).SM
    # Combine according to Eq.(treated_R0)
    return ((1 - pi_T) + phi_T * pi_T) * (1 - p)^2 * SM * R0_bar
end

"""
Compute the endemic equilibrium (EE) analytically by solving for I_H.
Returns NamedTuple of same form as compute_dfe, or nothing if R0 ≤ 1.
"""
function compute_endemic_equilibrium(par::Parameters)
    # R0 = compute_R0(par)
    # if R0 ≤ 1
    #     return nothing
    # end
    a = par.a
    b = par.b
    c = par.c
    m = par.m
    r = par.r
    g = par.g
    h1 = par.h1
    h2 = par.h2
    ε = par.ε
    p = par.p
    s1M = par.s1M
    s2M = par.s2M
    s1T = par.s1T
    s2T = par.s2T
    function g!(G, vars)
        SM, S2T = vars
        S1T = ε * p * a * (SM + S2T) / (h1 + g)
        G[1] = g + h2 * S2T - (((1 - p) * a * c * IH) + ε * p * a + g) * SM
        G[2] = h1 * S1T - ((ε * p * a * S2T + h2 + g) * S2T)
    end
    function f_IH(IH)
        function g!(G, vars)
            SM, S2T = vars
            S1T = ε * p * a * (SM + S2T) / (h1 + g)
            G[1] = g + h2 * S2T - (((1 - p) * a * c * IH) + ε * p * a + g) * SM
            G[2] = h1 * S1T - ((ε * p * a * S2T + h2 + g) * S2T)
        end
        sol = nlsolve(g!, [par.g, 0.0])
        SM, S2T = sol.zero
        E1M = (1 - p) * a * c * IH * SM / (ε * p * a + s1M + g)
        E2M = s1M * E1M / (s2M + g)
        IM = s2M * E2M / g
        E1T = ε * p * a * E1M / (s1T + g)
        E2T = s1T * E1T / (s2T + g)
        IT = s2T * E2T / g
        return m * (1 - p) * a * b * (IM + IT) * (1 - IH) - r * IH
    end
    lo, hi = 1e-6, 1 - 1e-6
    for _ in 1:50
        mid = (lo + hi) / 2
        if f_IH(mid) > 0
            hi = mid
        else
            lo = mid
        end
    end
    IH = (lo + hi) / 2
    SH = 1 - IH
    sol = nlsolve((G, v) -> g!(G, v), [par.g, 0.0])
    SM, S2T = sol.zero
    S1T = ε * p * a * (SM + S2T) / (h1 + g)
    E1M = (1 - p) * a * c * IH * SM / (ε * p * a + s1M + g)
    E2M = s1M * E1M / (s2M + g)
    IM = s2M * E2M / g
    E1T = ε * p * a * E1M / (s1T + g)
    E2T = s1T * E1T / (s2T + g)
    IT = s2T * E2T / g
    return (SH=SH, IH=IH, SM=SM, E1M=E1M, E2M=E2M, IM=IM,
        S1T=S1T, S2T=S2T, E1T=E1T, E2T=E2T, IT=IT)
end

"""
ODE system for numerical integration of the full model.
State vector: [SH,...,IT]
"""
function malaria_ode!(du, u, par::Parameters, t)
    SH, IH, SM, E1M, E2M, IM, S1T, S2T, E1T, E2T, IT = u
    a = par.a
    b = par.b
    c = par.c
    m = par.m
    r = par.r
    g = par.g
    h1 = par.h1
    h2 = par.h2
    ε = par.ε
    p = par.p
    s1M = par.s1M
    s2M = par.s2M
    s1T = par.s1T
    s2T = par.s2T
    du[1] = -m * (1 - p) * a * b * (IM + IT) * SH + r * IH
    du[2] = m * (1 - p) * a * b * (IM + IT) * SH - r * IH
    du[3] = g + h2 * S2T - (((1 - p) * a * c * IH) + ε * p * a + g) * SM
    du[4] = (1 - p) * a * c * IH * SM - (ε * p * a + s1M + g) * E1M
    du[5] = s1M * E1M - (s2M + g) * E2M
    du[6] = s2M * E2M - g * IM
    du[7] = ε * p * a * (SM + S2T) - (h1 + g) * S1T
    du[8] = h1 * S1T - ((ε * p * a) * S2T + h2 + g) * S2T
    du[9] = ε * p * a * E1M - (s1T + g) * E1T
    du[10] = s1T * E1T - (s2T + g) * E2T
    du[11] = s2T * E2T - g * IT
end

"""
Numerically compute the endemic equilibrium by integrating until convergence.
"""
function compute_ee_numerical(par::Parameters;
    u0=[0.99, 1e-3, 1.0, zeros(9)...],
    tspan=(0.0, 1.0),
    Nlast=100, tol=1e-6)
    buf = Vector{Vector{Float64}}()
    cur = copy(u0)
    for _ in 1:Nlast
        sol = solve(ODEProblem(malaria_ode!, cur, tspan, par), RK4())
        cur = sol.u[end]
        push!(buf, deepcopy(cur))
    end
    while true
        M = hcat(buf...)
        slopes = map(eachrow(M)) do row
            cov(collect(1:length(row)), row) / var(row)
        end
        if all(abs.(slopes) .< tol)
            break
        end
        sol = solve(ODEProblem(malaria_ode!, cur, tspan, par), RK4())
        cur = sol.u[end]
        popfirst!(buf)
        push!(buf, deepcopy(cur))
    end
    return vec(mean(hcat(buf...); dims=2))
end

end # module
