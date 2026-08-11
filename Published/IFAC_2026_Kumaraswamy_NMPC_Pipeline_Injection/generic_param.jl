##############################
# MODEL PARAMETERS  
##############################

# CO2 PHYSICAL PROPERTIES 
# RHO - kg/m3, Density
# EPS - m, Epsilon (Roughness)
# G - m/s2, Acceleration due to gravity
# MW - kg/mol,  Molecular Weight 
# R - m3.bar/(K.mol), Gas Constant 
# Tref - K, Reference Temperature 
# Tc - K, Critical Temperature 
# Pc - bar, Critical Pressure 

RHO = 927                                       
EPS = 20*(10^-6)                                
G = 9.81          
MW = 0.044 
R = 8.314*(10^-5)
Tref = 273.15 + 10 
Tc   = 304.1282           
Pc   = 73.773           

# CO2 PHYSICAL PROPERTIES VARIATION WITH PRESSURE - REGRESSION 

# Pressure: Pa 
# Viscosity: kg/(m.s)
# Compressibility Factor, Z: Dimensionless 
# Density: kg/m3
# Compressibility Values: 1/Pa 

P_vals = collect(50e5:1e5:300e5)
mu_vals = [PropsSI("viscosity", "T", Tref, "P", P, "CO2") for P in P_vals]
Z_vals = [PropsSI("Z", "T", Tref, "P", P, "CO2") for P in P_vals]
rho_vals = [PropsSI("D", "T", Tref, "P", P, "CO2") for P in P_vals]  # Density in kg/m³
compressibility_vals = [PropsSI("ISOTHERMAL_COMPRESSIBILITY", "T", Tref, "P", P, "CO2") for P in P_vals] # Compressibility in 1/ Pa

poly2(x, p) =  p[1] .+ p[2].*x .+ p[3].*x.^2 

# Fit the model using LsqFit.jl
fit_result_mu = curve_fit(poly2, P_vals, mu_vals, [1.0, 0.0, 0.0])
fit_result_rho = curve_fit(poly2, P_vals, rho_vals, [1.0, 0.0, 0.0])
fit_result_compressibility = curve_fit(poly2, P_vals, compressibility_vals, [1.0, 0.0, 0.0])
fit_result_Z = curve_fit(poly2, P_vals, Z_vals, [1.0, 0.0, 0.0])

# CAVITATION INDEX PARAMETERS 
co2_vapor_pressure = 50 # bar

##############################
# MODEL VARIABLES LIST   
##############################


model_vars = ["p_in", 
        "p_out",    
        "w_in",
        "w_out",
        "dpdz_in",
        "dpdz_out",
        "Z",
        "Z_in",
        "Z_out",
        "rho",
        "rho_in",
        "rho_out",
        "reynold_in",
        "reynold_out",
        "fdarcy_in",
        "fdarcy_out",
        "choke_vlv_op",
        "deltap_wellhead_choke",
        "p_node",
        "p_average",
        "rho_pump_inlet",
        "rho_pump_outlet",
        "pump_work",
        "cavitation_idx",
        "head_meters",
        "q_pump",
        "speed_pump",
        "injectivity",
        "p_reservoir"
]

