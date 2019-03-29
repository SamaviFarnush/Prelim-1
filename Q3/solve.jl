include("data.jl")
include("Flux.jl")

using PyPlot
using DelimitedFiles

pygui()

k_Ex = v_x/gene_length;                         #1/hr
k_El = v_l/peptide_length;                      #1/hr

number_of_simulation_points = 1000;
inducer_array = collect(exp10.(range(-4,stop=1,length=number_of_simulation_points)));
simulated_mRNA_concentration = zeros(length(inducer_array),4);
protein_flux = zeros(number_of_simulation_points);
flux = zeros(number_of_simulation_points,number_of_fluxes);

for step_index = 1:number_of_simulation_points

    inducer = inducer_array[step_index]
    f_binding = (inducer^n)/(K^n+(inducer)^n)
    u_value = (W1+W2*f_binding)/(1+W1+W2*f_binding)
    # compute the mRNA level -
    mRNA_level      = k_Ex*R_x*Gp/(K_x*tau_x+(tau_x+1)*Gp)*u_value;
    mRNA            = mRNA_level/kd_x;
    protein_level   = k_El*R_l*mRNA/(K_l*tau_l+(tau_l+1)*mRNA);

    objective_coefficient_array[5] = 1;

    default_bounds_array[7:end,2] .=  100000;
    default_bounds_array[7:end,1] .= -100000;
    default_bounds_array[1:6,2]   .=  100000;
    default_bounds_array[1:6,1]   .=  0;
    default_bounds_array[2,2]      =  mRNA_level;
    default_bounds_array[2,1]      =  mRNA_level;
    default_bounds_array[5,2]      =  protein_level;
    default_bounds_array[5,1]      =  0;

    (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag)=calculate_optimal_flux_distribution(stoichiometric_matrix, default_bounds_array, species_bounds_array, objective_coefficient_array)

    # cache -
    simulated_mRNA_concentration[step_index,1] = inducer
    simulated_mRNA_concentration[step_index,2] = mRNA_level
    simulated_mRNA_concentration[step_index,3] = u_value
    simulated_mRNA_concentration[step_index,4] = protein_level

    protein_production = objective_value;
    protein_flux[step_index] = protein_production;
    flux[step_index,:].=calculated_flux_array;
end
# 1d Make a plot - mRNA (y-axis) versus I (x-axis) -
semilogx(simulated_mRNA_concentration[:,1],protein_flux[1:end],color="blue",lw=1)

# label the axes -
xlabel("Inducer [mM]",fontsize=12)
ylabel("Optimized protein production [uM/hr]",fontsize=12)
title("Cell Free Protein Synthesis")
savefig("Q3.pdf")
gcf()
writedlm("flux.csv",flux,",")
