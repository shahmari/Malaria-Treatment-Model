using DifferentialEquations, LinearAlgebra, Plots, ProgressMeter, LaTeXStrings

include("basic_reproduction_number_generalized.jl")

using .BasicReproductionNumberGeneralized

# Define a function that returns a plausible parameter set with 2 latent stages.
function default_parameters_gen()
    return BasicReproductionNumberGeneralized.Parameters(
        0.2,    # a: Biting rate
        0.5,    # b: Human-to-mosquito transmission probability
        0.5,    # c: Mosquito-to-human transmission probability
        20.0,   # m: Mosquito-to-human ratio
        0.01,   # r: Human recovery rate (day⁻¹)
        0.12,   # g: Mosquito death rate (day⁻¹)
        0.1,    # h: Treatment waning rate (day⁻¹)
        0.1,    # t_rate: Treatment encounter rate (day⁻¹)
        2,      # N_LM: Number of latent stages for untreated mosquitoes
        2,      # N_LT: Number of latent stages for treated mosquitoes
        0.2,    # sigma_M: Progression rate for untreated latent stages (day⁻¹)
        0.1     # sigma_T: Progression rate for treated latent stages (day⁻¹)
    )
end

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
t_range_rev = LinRange(1e-10, 100.0, 500)
h_range_rev = LinRange(1e-10, 100.0, 500)

p.N_LM = 1
p.N_LT = 1

@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.h = 1 / h_range_rev[i]
    p.t_rate = 1 / t_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

txh_PLT = contour(t_range, h_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
txh_rev_PLT = contour(t_range_rev, h_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(txh_PLT, "../fig/gen_model/R0_rates_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(txh_PLT, "../fig/gen_model/R0_rates_txh_$(p.N_LM)x$(p.N_LT).pdf")
savefig(txh_rev_PLT, "../fig/gen_model/R0_periods_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(txh_rev_PLT, "../fig/gen_model/R0_periods_txh_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
t_range_rev = LinRange(1e-10, 100.0, 500)
h_range_rev = LinRange(1e-10, 100.0, 500)

p.N_LM = 2
p.N_LT = 2

@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.h = 1 / h_range_rev[i]
    p.t_rate = 1 / t_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

txh_PLT = contour(t_range, h_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
txh_rev_PLT = contour(t_range_rev, h_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(txh_PLT, "../fig/gen_model/R0_rates_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(txh_PLT, "../fig/gen_model/R0_rates_txh_$(p.N_LM)x$(p.N_LT).pdf")
savefig(txh_rev_PLT, "../fig/gen_model/R0_periods_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(txh_rev_PLT, "../fig/gen_model/R0_periods_txh_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
t_range_rev = LinRange(1e-10, 100.0, 500)
h_range_rev = LinRange(1e-10, 100.0, 500)

p.N_LM = 5
p.N_LT = 5

@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.h = 1 / h_range_rev[i]
    p.t_rate = 1 / t_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

txh_PLT = contour(t_range, h_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
txh_rev_PLT = contour(t_range_rev, h_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(txh_PLT, "../fig/gen_model/R0_rates_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(txh_PLT, "../fig/gen_model/R0_rates_txh_$(p.N_LM)x$(p.N_LT).pdf")
savefig(txh_rev_PLT, "../fig/gen_model/R0_periods_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(txh_rev_PLT, "../fig/gen_model/R0_periods_txh_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
t_range_rev = LinRange(1e-10, 100.0, 500)
h_range_rev = LinRange(1e-10, 100.0, 500)

p.N_LM = 10
p.N_LT = 10

@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.h = 1 / h_range_rev[i]
    p.t_rate = 1 / t_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

txh_PLT = contour(t_range, h_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
txh_rev_PLT = contour(t_range_rev, h_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(txh_PLT, "../fig/gen_model/R0_rates_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(txh_PLT, "../fig/gen_model/R0_rates_txh_$(p.N_LM)x$(p.N_LT).pdf")
savefig(txh_rev_PLT, "../fig/gen_model/R0_periods_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(txh_rev_PLT, "../fig/gen_model/R0_periods_txh_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
SM_range = LinRange(0.0, 3.0, 500)
ST_range = LinRange(0.0, 3.0, 500)
SM_range_rev = LinRange(1e-10, 10.0, 500)
ST_range_rev = LinRange(1e-10, 10.0, 500)

p.N_LM = 1
p.N_LT = 1

@showprogress for i in 1:500, j in 1:500
    p.sigma_T = ST_range[i]
    p.sigma_M = SM_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.sigma_T = 1 / ST_range_rev[i]
    p.sigma_M = 1 / SM_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

SMxST_PLT = contour(SM_range, ST_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
SMxST_rev_PLT = contour(SM_range_rev, ST_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency period, " * L"S_M (day)", ylabel="Treated latency period, " * L"S_T (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).pdf")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
SM_range = LinRange(0.0, 3.0, 500)
ST_range = LinRange(0.0, 3.0, 500)
SM_range_rev = LinRange(1e-10, 10.0, 500)
ST_range_rev = LinRange(1e-10, 10.0, 500)

p.N_LM = 2
p.N_LT = 2

@showprogress for i in 1:500, j in 1:500
    p.sigma_T = ST_range[i]
    p.sigma_M = SM_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.sigma_T = 1 / ST_range_rev[i]
    p.sigma_M = 1 / SM_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

SMxST_PLT = contour(SM_range, ST_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
SMxST_rev_PLT = contour(SM_range_rev, ST_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency period, " * L"S_M (day)", ylabel="Treated latency period, " * L"S_T (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).pdf")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
SM_range = LinRange(0.0, 3.0, 500)
ST_range = LinRange(0.0, 3.0, 500)
SM_range_rev = LinRange(1e-10, 10.0, 500)
ST_range_rev = LinRange(1e-10, 10.0, 500)

p.N_LM = 5
p.N_LT = 5

@showprogress for i in 1:500, j in 1:500
    p.sigma_T = ST_range[i]
    p.sigma_M = SM_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.sigma_T = 1 / ST_range_rev[i]
    p.sigma_M = 1 / SM_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

SMxST_PLT = contour(SM_range, ST_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
SMxST_rev_PLT = contour(SM_range_rev, ST_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency period, " * L"S_M (day)", ylabel="Treated latency period, " * L"S_T (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).pdf")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
SM_range = LinRange(0.0, 3.0, 500)
ST_range = LinRange(0.0, 3.0, 500)
SM_range_rev = LinRange(1e-10, 10.0, 500)
ST_range_rev = LinRange(1e-10, 10.0, 500)

p.N_LM = 10
p.N_LT = 10

@showprogress for i in 1:500, j in 1:500
    p.sigma_T = ST_range[i]
    p.sigma_M = SM_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.sigma_T = 1 / ST_range_rev[i]
    p.sigma_M = 1 / SM_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

SMxST_PLT = contour(SM_range, ST_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
SMxST_rev_PLT = contour(SM_range_rev, ST_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency period, " * L"S_M (day)", ylabel="Treated latency period, " * L"S_T (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).pdf")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
R0_Matrix_Rates = zeros(500, 500)
R0_Matrix_Periods = zeros(500, 500)
SM_range = LinRange(0.0, 3.0, 500)
ST_range = LinRange(0.0, 3.0, 500)
SM_range_rev = LinRange(1e-10, 10.0, 500)
ST_range_rev = LinRange(1e-10, 10.0, 500)

p.N_LM = 50
p.N_LT = 50

@showprogress for i in 1:500, j in 1:500
    p.sigma_T = ST_range[i]
    p.sigma_M = SM_range[j]
    R0_Matrix_Rates[i, j] = compute_R0_generalized(p)
    p.sigma_T = 1 / ST_range_rev[i]
    p.sigma_M = 1 / SM_range_rev[j]
    R0_Matrix_Periods[i, j] = compute_R0_generalized(p)
end

SMxST_PLT = contour(SM_range, ST_range, R0_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency rate, " * L"S_M (day^{-1})", ylabel="Treated latency rate, " * L"S_T (day^{-1})", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))
SMxST_rev_PLT = contour(SM_range_rev, ST_range_rev, R0_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Untreated latency period, " * L"S_M (day)", ylabel="Treated latency period, " * L"S_T (day)", title="Basic Reproduction Number, " * L"R_0 \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.0, 13.0), size=(600, 500))

savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_PLT, "../fig/gen_model/R0_rates_SMxST_$(p.N_LM)x$(p.N_LT).pdf")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).png")
savefig(SMxST_rev_PLT, "../fig/gen_model/R0_periods_SMxST_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
IH_Matrix_Rates = zeros(500, 500)
IH_Matrix_Periods = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
t_range_rev = LinRange(1e-10, 100.0, 500)
h_range_rev = LinRange(1e-10, 100.0, 500)

p.N_LM = 1
p.N_LT = 1

@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    EE_Solution = compute_endemic_equilibrium_generalized(p)
    IH_Matrix_Rates[i, j] = EE_Solution[2]
    p.h = 1 / h_range_rev[i]
    p.t_rate = 1 / t_range_rev[j]
    EE_Solution = compute_endemic_equilibrium_generalized(p)
    IH_Matrix_Periods[i, j] = EE_Solution[2]
end

IH_PLT = contour(t_range, h_range, IH_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Infectious Humans, " * L"I_H \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.5, 1.0), size=(600, 500))
IH_rev_PLT = contour(t_range_rev, h_range_rev, IH_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Infectious Humans, " * L"I_H \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.5, 1.0), size=(600, 500))

savefig(IH_PLT, "../fig/gen_model/IH_rates_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(IH_PLT, "../fig/gen_model/IH_rates_txh_$(p.N_LM)x$(p.N_LT).pdf")
savefig(IH_rev_PLT, "../fig/gen_model/IH_periods_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(IH_rev_PLT, "../fig/gen_model/IH_periods_txh_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
IH_Matrix_Rates = zeros(500, 500)
IH_Matrix_Periods = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
t_range_rev = LinRange(1e-10, 100.0, 500)
h_range_rev = LinRange(1e-10, 100.0, 500)

p.N_LM = 2
p.N_LT = 2

@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    EE_Solution = compute_endemic_equilibrium_generalized(p)
    IH_Matrix_Rates[i, j] = EE_Solution[2]
    p.h = 1 / h_range_rev[i]
    p.t_rate = 1 / t_range_rev[j]
    EE_Solution = compute_endemic_equilibrium_generalized(p)
    IH_Matrix_Periods[i, j] = EE_Solution[2]
end

IH_PLT = contour(t_range, h_range, IH_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Infectious Humans, " * L"I_H \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.5, 1.0), size=(600, 500))
IH_rev_PLT = contour(t_range_rev, h_range_rev, IH_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Infectious Humans, " * L"I_H \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.5, 1.0), size=(600, 500))

savefig(IH_PLT, "../fig/gen_model/IH_rates_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(IH_PLT, "../fig/gen_model/IH_rates_txh_$(p.N_LM)x$(p.N_LT).pdf")
savefig(IH_rev_PLT, "../fig/gen_model/IH_periods_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(IH_rev_PLT, "../fig/gen_model/IH_periods_txh_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()
IH_Matrix_Rates = zeros(500, 500)
IH_Matrix_Periods = zeros(500, 500)
t_range = LinRange(0.0, 1.0, 500)
h_range = LinRange(0.0, 1.0, 500)
t_range_rev = LinRange(1e-10, 100.0, 500)
h_range_rev = LinRange(1e-10, 100.0, 500)

p.N_LM = 5
p.N_LT = 5

@showprogress for i in 1:500, j in 1:500
    p.h = h_range[i]
    p.t_rate = t_range[j]
    EE_Solution = compute_endemic_equilibrium_generalized(p; tol=1e-5, max_iter=100000)
    IH_Matrix_Rates[i, j] = EE_Solution[2]
    p.h = 1 / h_range_rev[i]
    p.t_rate = 1 / t_range_rev[j]
    EE_Solution = compute_endemic_equilibrium_generalized(p; tol=1e-5, max_iter=100000)
    IH_Matrix_Periods[i, j] = EE_Solution[2]
end

IH_PLT = contour(t_range, h_range, IH_Matrix_Rates, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter rate, t (day⁻¹)", ylabel="Wane rate, h (day⁻¹)", title="Infectious Humans, " * L"I_H \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.5, 1.0), size=(600, 500))
IH_rev_PLT = contour(t_range_rev, h_range_rev, IH_Matrix_Periods, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Treatment encounter period, t⁻¹ (day)", ylabel="Wane period, h⁻¹ (day)", title="Infectious Humans, " * L"I_H \ (N_{L,M} = %$(p.N_LM), N_{L,T} = %$(p.N_LT))", clims=(0.5, 1.0), size=(600, 500))

savefig(IH_PLT, "../fig/gen_model/IH_rates_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(IH_PLT, "../fig/gen_model/IH_rates_txh_$(p.N_LM)x$(p.N_LT).pdf")
savefig(IH_rev_PLT, "../fig/gen_model/IH_periods_txh_$(p.N_LM)x$(p.N_LT).png")
savefig(IH_rev_PLT, "../fig/gen_model/IH_periods_txh_$(p.N_LM)x$(p.N_LT).pdf")

p = default_parameters_gen()

N_LM_range = 1:100
N_LT_range = 1:100
R0_Matrix = zeros(length(N_LM_range), length(N_LT_range))

@showprogress for i in 1:100, j in 1:100
    p.N_LT = N_LT_range[i]
    p.N_LM = N_LM_range[j]
    R0_Matrix[i, j] = compute_R0_generalized(p)
end

R0_PLT = contour(N_LM_range, N_LT_range, R0_Matrix, c=cgrad(:thermal, rev=true), lc=:black, fill=true,
    xlabel="Number of untreated latent stages, " * L"N_{L,M}", ylabel="Number of treated latent stages, " * L"N_{L,T}", title="Basic Reproduction Number, " * L"R_0", clims=(0.0, 13.0), size=(600, 500))

savefig(R0_PLT, "../fig/gen_model/R0_NLMxNLT.png")
savefig(R0_PLT, "../fig/gen_model/R0_NLMxNLT.pdf")