using DifferentialEquations, LinearAlgebra, Plots, Plots.Measures, ProgressMeter, LaTeXStrings

Plots.default(titlefontsize=16, tickfontsize=14, labelfontsize=16, legendfontsize=10, rightmargin=5mm,
    fontfamily="Sans-serif", frame=:box, label=nothing)

cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file

include("basic_reproduction_number_new_model.jl")

using .MalariaTreatmentModel

# 1) Set up a default parameter constructor
function default_parameters()
    return Parameters(
        0.028,  # a: biting rate
        0.50,   # b: human→mosquito trans. prob.
        0.50,   # c: mosquito→human trans. prob.
        20.0,   # m: mosquito→human ratio
        0.01,   # r: human recovery rate
        0.12,   # g: mosquito death rate
        0.10,   # h1: treated S1→S2 waning rate
        0.10,   # h2: treated S2→SM waning rate
        1.0,    # ε: treatment efficacy factor
        1.0,    # p: treatment coverage
        0.20,   # s1M: untreated latent progress. rate (stage 1→2)
        0.20,   # s2M: untreated latent progress. rate (stage 2→I)
        0.10,   # s1T:   treated latent progress. rate (stage 1→2)
        0.10,   # s2T:   treated latent progress. rate (stage 2→I)
        2       # L: number of latent stages
    )
end

# 2) Define our sweep grid
N = 500
ε_range = LinRange(0.0, 1.0, N)
p_range = LinRange(0.0, 1.0, N)

# 3) Pre-allocate storage
R0_mat = zeros(N, N)
psi_mat = zeros(N, N)
SM_mat = zeros(N, N)
cov2_mat = zeros(N, N)
mix_mat = zeros(N, N)

# 4) Loop and compute
p = default_parameters()
@showprogress for i in 1:N, j in 1:N
    p.ε = ε_range[i]
    p.p = p_range[j]
    # closed-form results
    vals = compute_R0_closed_form(p)
    R0, R0_bar, pi_T, phi_T = vals.R0, vals.R0_bar, vals.pi_T, vals.phi_T

    # compute derived quantities
    ψ = R0 / R0_bar
    SM = compute_dfe_closed_form(p).SM
    c2 = (1 - p.p)^2
    mix = 1 - pi_T + phi_T * pi_T

    # store
    R0_mat[i, j] = R0
    psi_mat[i, j] = ψ
    SM_mat[i, j] = SM
    cov2_mat[i, j] = c2
    mix_mat[i, j] = mix
end

# 5) Plot heatmaps
plots = Dict{String,Any}()

plots["R0"] = contour(
    ε_range, p_range, R0_mat',
    c=cgrad(:thermal, rev=true), fill=true, lc=:black,
    xlabel=L"\varepsilon", ylabel=L"p",
    title=L"R_0(\varepsilon,p)",
    size=(750, 600),
    rightmargin=10mm
)

plots["psi"] = contour(
    ε_range, p_range, psi_mat',
    c=cgrad(:thermal, rev=true), fill=true, lc=:black,
    xlabel=L"\varepsilon", ylabel=L"p",
    title=L"\psi = R_0 / R_{0,\mathrm{bar}}",
    size=(750, 600),
    rightmargin=10mm
)

plots["SMstar"] = contour(
    ε_range, p_range, SM_mat',
    c=cgrad(:thermal, rev=true), fill=true, lc=:black,
    xlabel=L"\varepsilon", ylabel=L"p",
    title=L"S_M^*(\varepsilon,p)",
    size=(750, 600),
    rightmargin=10mm
)

plots["coverage"] = contour(
    ε_range, p_range, cov2_mat',
    c=cgrad(:thermal, rev=true), fill=true, lc=:black,
    xlabel=L"\varepsilon", ylabel=L"p",
    title=L"(1 - p)^2",
    size=(750, 600),
    rightmargin=10mm
)

plots["mix"] = contour(
    ε_range, p_range, mix_mat',
    c=cgrad(:thermal, rev=true), fill=true, lc=:black,
    xlabel=L"\varepsilon", ylabel=L"p",
    title=L"1 - \pi_T + \phi_T\,\pi_T",
    size=(750, 600),
    rightmargin=10mm
)

# 6) Save them
for (name, plt) in plots
    savefig(plt, "../fig/nm_heatmap_$(name).png")
    savefig(plt, "../fig/nm_heatmap_$(name).pdf")
end
