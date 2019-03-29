Gp              = 5*10^-3;       #uM
cell_volume     = 15*10^-6;      #L
R_x             = 0.15;          #uM
R_l             = 1.6;           #uM
v_x             = 60*3600;       #nt/h
v_l             = 16.5*3600;     #aa/h
K_x             = 0.3;           #uM
K_l             = 57;            #uM
tau_x           = 2.7;           #dimensionless
tau_l           = 0.8;           #dimensionless
kd_x            = 8.35;          #1/h
kd_l            = 9.9*10^-3;     #1/h
peptide_length  = 308;           #aa
gene_length     = 924;           #nt
W1              = 0.26;
W2              = 300;
K               = 0.3;           #mM
n               = 1.5;

a=1;
n=1;

stoichiometric_matrix = [
#    v1  v2  v3  v4  v5  v6  b1  b2  b3  b4  b5  b6  b7  b8  b9
   -1.0   1   0   0   0   0   0   0   0   0   0   0   0   0   0     #Gp
     -1   1   0   0   0   0   0   0   0   0   0   0   0   0   0     #RNAP
      1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0     #G*
      0  -n   0   0   0   0   0   1   0   0   0   0   0   0   0     #NTP
      0   1  -1  -1   1   0   0   0   0   0   0   0   0   0   0     #mRNA
      0  2n   0   0  2a   2   0   0   0   0   0   0   0   0  -1     #Pi
      0   0   1   0   0   0   0   0   0  -1   0   0   0   0   0     #NMP
      0   0   0  -1   1   0   0   0   0   0   0   0   0   0   0     #rib
      0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0     #rib*
      0   0   0   0  -a   1   0   0   0   0   0   0   0   0   0     #AAtRNA
      0   0   0   0 -2a   0   0   0   0   0   0   0   1   0   0     #GTP
      0   0   0   0   a  -1   0   0   0   0   0   0   0   0   0     #tRNA
      0   0   0   0  2a   0   0   0   0   0   0   0   0  -1   0     #GDP
      0   0   0   0   1   0   0   0  -1   0   0   0   0   0   0     #protein
      0   0   0   0   0  -1   1   0   0   0   0   0   0   0   0     #AA
      0   0   0   0   0  -1   0   0   0   0   1   0   0   0   0     #ATP
      0   0   0   0   0   1   0   0   0   0   0  -1   0   0   0     #AMP
]

(number_of_species,number_of_fluxes) = size(stoichiometric_matrix);

default_bounds_array = zeros(number_of_fluxes,2);

## maximization function (in this case which is urea production flux b4)
objective_coefficient_array = zeros(number_of_fluxes);

##RHS of S*v==0
species_bounds_array = zeros(number_of_species,2);
