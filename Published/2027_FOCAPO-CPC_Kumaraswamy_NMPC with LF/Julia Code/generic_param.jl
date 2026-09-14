
##############################
# MODEL PARAMETERS
##############################
# Parameters shared by every case. Case-specific settings live in case_<name>.jl, which must be
# included before this file because the property fits below are evaluated at include time.

# CO2 PHYSICAL PROPERTIES
# RHO - kg/m3, nominal density (reference value only; the models use the pressure-dependent fit below)
# EPS - m, pipe wall roughness
# MU - kg/m.s, nominal dynamic viscosity (reference value only; see the fit below)
# G - m/s2, Acceleration due to gravity
# MW - kg/mol,  Molecular Weight
# R - m3.bar/(K.mol), Gas Constant
# Tref - K, Reference Temperature. The models are isothermal, so every property is evaluated here
# Tc - K, Critical Temperature
# Pc - bar, Critical Pressure

RHO = 927                                       
EPS = 20*(10^-6)                                
MU = 13.7*(10^-6)                               
G = 9.81          
MW = 0.044 
R = 8.314*(10^-5)
Tref = 273.15 + 10 
Tc   = 304.1282           
Pc   = 73.773           

# CO2 PHYSICAL PROPERTIES VARIATION WITH PRESSURE - REGRESSION
# CoolProp is sampled on a 50-300 bar grid at Tref, then each property is fitted with a quadratic
# in pressure. The NMPC calls the fits rather than CoolProp because the fits are smooth and
# differentiable, which is what Ipopt needs.
#
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

# Quadratic in pressure; `p` is the coefficient vector returned by curve_fit as fit_result_*.param
poly2(x, p) =  p[1] .+ p[2].*x .+ p[3].*x.^2

# Fit the model using LsqFit.jl. Only fit_result_mu and fit_result_Z are used by the NMPC and the
# steady-state model; density is recovered from Z through the real gas law rather than from
# fit_result_rho, and fit_result_compressibility is kept for reference.
fit_result_mu = curve_fit(poly2, P_vals, mu_vals, [1.0, 0.0, 0.0])
fit_result_rho = curve_fit(poly2, P_vals, rho_vals, [1.0, 0.0, 0.0])
fit_result_compressibility = curve_fit(poly2, P_vals, compressibility_vals, [1.0, 0.0, 0.0])
fit_result_Z = curve_fit(poly2, P_vals, Z_vals, [1.0, 0.0, 0.0])

# CAVITATION INDEX PARAMETERS
co2_vapor_pressure = 50 # bar, vapour pressure of CO2 at Tref, used as the cavitation reference

# COLLOCATION PARAMETERS
# The dynamic model is discretised in time by orthogonal collocation on finite elements, one
# element per sample of length DT. D is the number of collocation points per element: d = 1 is
# the element start (carried over from the previous sample) and d = 2:D are the interior points
# Adot is the (D x D-1) matrix of Lagrange basis derivatives at the collocation points, used to
# turn the mass balance ODE into the algebraic collocation equations.
D = 4
Adot =

  [-4.139387691339813   1.739387691339811  -3.000000000000002
   3.224744871391587  -3.567840084690404   5.531972647421805
   1.167840084690405   0.775255128608409  -7.531972647421807
  -0.253197264742181   1.053197264742181   5.000000000000000]

##############################
# MODEL VARIABLES LIST
##############################
# Variables that exist in both model_ss.jl and model_mpc.jl under the same name. 
# The list is used for three things: 
# 1. Extracts an initial guess from the SS model for the first NMPC solve (case notebook)
# 2. Saves the full-horizon predictions
# 3. Shifts the solution forward one sample to warm-start the next solve (run_mpc.jl). 
# Note: A variable added to either model must be added here to be carried.

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
        "p_reservoir",
]

# Bias-corrected outputs. These exist only in the NMPC, so they are saved with the predictions but
# are not part of the warm start and have no steady-state counterpart.
corrected_terms =
[        "w_in_corrected",
        "w_in_corrected_mainline",
        "p_node_corrected"]