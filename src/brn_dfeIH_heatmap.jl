using Plots, LaTeXStrings, Plots.Measures, ProgressMeter

Plots.default(titlefontsize=16, tickfontsize=14, labelfontsize=16, legendfontsize=10, rightmargin=5mm,
    fontfamily="Sans-serif", frame=:box, label=nothing)

include("basic_reproduction_number.jl")

using .BasicReproductionNumber

FigureDir = "../fig/" # Directory for storing the figures
cd(dirname(@__FILE__)) # Changing the current directory to the directory of the file
println("Directories and modules have been set...\n")

# Define a function that returns a plausible parameter set
function default_parameters()
    return Parameters(
        0.2,    # Biting rate
        0.5,    # Human-to-mosquito transmission
        0.5,    # Mosquito-to-human transmission
        20.0,   # Mosquito-to-human ratio
        0.01,    # Human recovery rate
        0.12,   # Mosquito death rate
        0.1,   # Treatment waning rate
        0.1, # Treatment encounter rate
        0.2,  # Untreated E1 -> E2 rate
        0.2,  # Untreated E2 -> I rate
        0.1, # Treated E1 -> E2 rate
        0.1   # Treated E2 -> I rate
    )
end

# Calculate the basic reproduction number
p = default_parameters()
R0_Matrix_m = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    R0_Matrix_m[i, j] = compute_R0(p)
end

p = default_parameters()
R0_Matrix_rev_m = zeros(500, 500)
rev_t_range = LinRange(0.001, 100.0, 500)
rev_h_range = LinRange(0.001, 100.0, 500)
@showprogress for i in 1:500, j in 1:500
    p.h = 1 / rev_h_range[i]
    p.t_rate = 1 / rev_t_range[j]
    R0_Matrix_rev_m[i, j] = compute_R0(p)
end

p = default_parameters()
R0_Matrix_s = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    R0_Matrix_s[i, j] = compute_R0_symbolic(p)
end

p = default_parameters()
R0_Matrix_rev_s = zeros(500, 500)
rev_t_range = LinRange(0.001, 100.0, 500)
rev_h_range = LinRange(0.001, 100.0, 500)
@showprogress for i in 1:500, j in 1:500
    p.h = 1 / rev_h_range[i]
    p.t_rate = 1 / rev_t_range[j]
    R0_Matrix_rev_s[i, j] = compute_R0_symbolic(p)
end

# Plot the heatmap
R0_Matrix_m_PLT = contour(t_range, h_range, R0_Matrix_m, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Basic Reproduction Number, R₀", size=(600, 500))

R0_Matrix_rev_m_PLT = contour(rev_t_range, rev_h_range, R0_Matrix_rev_m, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Basic Reproduction Number, R₀", size=(600, 500))

R0_Matrix_s_PLT = contour(t_range, h_range, R0_Matrix_s, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Basic Reproduction Number, R₀", size=(600, 500))

R0_Matrix_rev_s_PLT = contour(rev_t_range, rev_h_range, R0_Matrix_rev_s, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Basic Reproduction Number, R₀", size=(600, 500))

# Save the figure
savefig(R0_Matrix_m_PLT, FigureDir * "brn_txh_heatmap_m.png")
savefig(R0_Matrix_m_PLT, FigureDir * "brn_txh_heatmap_m.pdf")

savefig(R0_Matrix_rev_m_PLT, FigureDir * "brn_txh_heatmap_rev_m.png")
savefig(R0_Matrix_rev_m_PLT, FigureDir * "brn_txh_heatmap_rev_m.pdf")

savefig(R0_Matrix_s_PLT, FigureDir * "brn_txh_heatmap_s.png")
savefig(R0_Matrix_s_PLT, FigureDir * "brn_txh_heatmap_s.pdf")

savefig(R0_Matrix_rev_s_PLT, FigureDir * "brn_txh_heatmap_rev_s.png")
savefig(R0_Matrix_rev_s_PLT, FigureDir * "brn_txh_heatmap_rev_s.pdf")


# Heatmap for other parameters

p = default_parameters()
R0_Matrix = zeros(500, 500)
S_M_range = LinRange(0.0, 2.0, 500)
S_T_range = LinRange(0.0, 2.0, 500)
@showprogress for i in 1:500, j in 1:500
    p.s1T = S_T_range[i]
    p.s2T = S_T_range[i]
    p.s1M = S_M_range[j]
    p.s2M = S_M_range[j]
    R0_Matrix[i, j] = compute_R0_symbolic(p)
end

R0_PLT = contour(S_M_range, S_T_range, R0_Matrix, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", title="Basic Reproduction Number, R₀", size=(600, 500))
savefig(R0_PLT, "../fig/R0_STxSM.pdf")
savefig(R0_PLT, "../fig/R0_STxSM.png")

p = default_parameters()
R0_Matrix_rev = zeros(500, 500)
S_M_range_rev = LinRange(0.00, 10.0, 500)
S_T_range_rev = LinRange(0.00, 10.0, 500)
@showprogress for i in 1:500, j in 1:500
    p.s1T = 1 / S_T_range_rev[i]
    p.s2T = 1 / S_T_range_rev[i]
    p.s1M = 1 / S_M_range_rev[j]
    p.s2M = 1 / S_M_range_rev[j]
    R0_Matrix_rev[i, j] = compute_R0_symbolic(p)
end

R0_PLT = contour(S_M_range_rev, S_T_range_rev, R0_Matrix_rev, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency period, " * L"S_M (day)", ylabel="Treated latency period, " * L"S_T (day)", title="Basic Reproduction Number, R₀", size=(600, 500))
savefig(R0_PLT, "../fig/R0_STxSM_rev.pdf")
savefig(R0_PLT, "../fig/R0_STxSM_rev.png")

p = default_parameters()
R0_Matrix = zeros(500, 500)
S_M_range = LinRange(0.0, 0.1, 500)
S_T_range = LinRange(0.0, 0.1, 500)
@showprogress for i in 1:500, j in 1:500
    p.s1T = S_T_range[i]
    p.s2T = S_T_range[i]
    p.s1M = S_M_range[j]
    p.s2M = S_M_range[j]
    R0_Matrix[i, j] = compute_R0_symbolic(p)
end

R0_PLT = contour(S_M_range, S_T_range, R0_Matrix, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", title="Basic Reproduction Number, R₀", size=(600, 500))
R0_Cut_Plt = contour(S_M_range, S_T_range, R0_Matrix .> 1.0, c=cgrad([:white, :green], rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", legend=false, size=(600, 500))

savefig(R0_PLT, "../fig/R0_STxSM_zoom.pdf")
savefig(R0_PLT, "../fig/R0_STxSM_zoom.png")

savefig(R0_Cut_Plt, "../fig/R0_STxSM_zoom_cut.pdf")
savefig(R0_Cut_Plt, "../fig/R0_STxSM_zoom_cut.png")

p = default_parameters()
Ih_Matrix = zeros(500, 500)
t_range = LinRange(0.0, 5.0, 500)
h_range = LinRange(0.0, 5.0, 500)
@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    EE_Solution = compute_endemic_equilibrium(p)
    I_h = EE_Solution[2]
    Ih_Matrix[i, j] = I_h
end

p = default_parameters()
Ih_Matrix_rev = zeros(500, 500)
rev_t_range = LinRange(0.0, 2.0, 500)
rev_h_range = LinRange(0.0, 2.0, 500)
@showprogress for i in 1:500, j in 1:500
    p.h = 1 / rev_h_range[i]
    p.t_rate = 1 / rev_t_range[j]
    EE_Solution = compute_endemic_equilibrium(p)
    I_h = EE_Solution[2]
    Ih_Matrix_rev[i, j] = I_h
end

Ih_PLT = contour(t_range, h_range, Ih_Matrix, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Endemic human infected frequency, " * L"I_H", size=(600, 500))

Ih_PLT_rev = contour(rev_t_range, rev_h_range, Ih_Matrix_rev, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Endemic human infected frequency, " * L"I_H", size=(600, 500))

savefig(Ih_PLT, "../fig/Ih_txh.pdf")
savefig(Ih_PLT, "../fig/Ih_txh.png")

savefig(Ih_PLT_rev, "../fig/Ih_txh_rev.pdf")
savefig(Ih_PLT_rev, "../fig/Ih_txh_rev.png")

p = default_parameters()
Ih_Matrix = zeros(500, 500)
S_M_range = LinRange(0.0, 0.2, 500)
S_T_range = LinRange(0.0, 0.2, 500)
@showprogress for i in 1:500, j in 1:500
    p.s1T = S_T_range[i]
    p.s2T = S_T_range[i]
    p.s1M = S_M_range[j]
    p.s2M = S_M_range[j]
    EE_Solution = compute_endemic_equilibrium(p)
    I_h = EE_Solution[2]
    Ih_Matrix[i, j] = I_h
end

Ih_PLT = contour(S_M_range, S_T_range, Ih_Matrix, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", title="Endemic human infected frequency, " * L"I_H", size=(600, 500))

savefig(Ih_PLT, "../fig/Ih_STxSM.pdf")
savefig(Ih_PLT, "../fig/Ih_STxSM.png")

p = default_parameters()
Ih_Matrix_rev = zeros(500, 500)
S_M_range_rev = LinRange(0.0, 100.0, 500)
S_T_range_rev = LinRange(0.0, 100.0, 500)

@showprogress for i in 1:500, j in 1:500
    p.s1T = 1 / S_T_range_rev[i]
    p.s2T = 1 / S_T_range_rev[i]
    p.s1M = 1 / S_M_range_rev[j]
    p.s2M = 1 / S_M_range_rev[j]
    EE_Solution = compute_endemic_equilibrium(p)
    I_h = EE_Solution[2]
    Ih_Matrix_rev[i, j] = I_h
end

Ih_PLT_rev = contour(S_M_range_rev, S_T_range_rev, Ih_Matrix_rev, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency period, " * L"S_M (day)", ylabel="Treated latency period, " * L"S_T (day)", title="Endemic human infected frequency, " * L"I_H", size=(600, 500))

savefig(Ih_PLT_rev, "../fig/Ih_STxSM_rev.pdf")
savefig(Ih_PLT_rev, "../fig/Ih_STxSM_rev.png")

println("Basic reproduction number heatmap has been plotted and saved...\n")