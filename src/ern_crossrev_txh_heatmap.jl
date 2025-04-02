using Plots, LaTeXStrings, Plots.Measures, ProgressMeter, JLD

Plots.default(titlefontsize=16, tickfontsize=14, labelfontsize=16, legendfontsize=10, rightmargin=20mm,
    fontfamily="Sans-serif", frame=:box, label=nothing)

include("effective_reproduction_number.jl")

using .EffectiveReproductionNumber

FigureDir = "../fig/" # Directory for storing the figures
DataDir = "../data/" # Directory for storing the data
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

# Calculate the effective reproduction number
try
    Re_Data = load(DataDir * "ern_matrix.jld")
    global Re_Matrix = Re_Data["Re_Matrix"]
    global t_range = Re_Data["t_range"]
    global h_range = Re_Data["h_range"]
    println("Effective reproduction number heatmap has been loaded...\n")
catch
    p = default_parameters()
    Matrix_Dim = 500
    global Re_Matrix = zeros(Matrix_Dim, Matrix_Dim)
    global t_range = LinRange(0.01, 1.0, Matrix_Dim)
    global h_range = LinRange(0.01, 1.0, Matrix_Dim)
    @showprogress for i in 1:Matrix_Dim, j in 1:Matrix_Dim
        p.h = h_range[i]
        p.t_rate = t_range[j]
        Re_Matrix[i, j] = compute_Re(p)
    end
    save(DataDir * "ern_matrix.jld", "Re_Matrix", Re_Matrix, "t_range", t_range, "h_range", h_range)
    println("Effective reproduction number matrix has been calculated and saved...\n")
end

try
    Re_Data_rev = load(DataDir * "ern_matrix_rev.jld")
    global Re_Matrix_rev = Re_Data_rev["Re_Matrix_rev"]
    global rev_t_range = Re_Data_rev["rev_t_range"]
    global rev_h_range = Re_Data_rev["rev_h_range"]
    println("Effective reproduction number heatmap (reverse) has been loaded...\n")
catch
    p = default_parameters()
    Matrix_Dim = 500
    global Re_Matrix_rev = zeros(Matrix_Dim, Matrix_Dim)
    global rev_t_range = LinRange(0.01, 20.0, Matrix_Dim)
    global rev_h_range = LinRange(0.01, 20.0, Matrix_Dim)
    @showprogress for i in 1:Matrix_Dim, j in 1:Matrix_Dim
        p.h = 1 / rev_h_range[i]
        p.t_rate = 1 / rev_t_range[j]
        Re_Matrix_rev[i, j] = compute_Re(p)
    end
    save(DataDir * "ern_matrix_rev.jld", "Re_Matrix_rev", Re_Matrix_rev, "rev_t_range", rev_t_range, "rev_h_range", rev_h_range)
    println("Effective reproduction number matrix (reverse) has been calculated and saved...\n")
end

# Plot the heatmap
Re_Matrix_PLT = heatmap(t_range, h_range, Re_Matrix, c=:thermal,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Effective Reproduction Number, R₀", size=(650, 500))

Re_Matrix_rev_PLT = heatmap(rev_t_range, rev_h_range, Re_Matrix_rev, c=:thermal,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Effective Reproduction Number, R₀", size=(650, 500))

# Save the figure
savefig(Re_Matrix_PLT, FigureDir * "ern_txh_heatmap.png")
savefig(Re_Matrix_PLT, FigureDir * "ern_txh_heatmap.pdf")

savefig(Re_Matrix_rev_PLT, FigureDir * "ern_txh_heatmap_rev.png")
savefig(Re_Matrix_rev_PLT, FigureDir * "ern_txh_heatmap_rev.pdf")
println("Effective reproduction number heatmap has been plotted and saved...\n")