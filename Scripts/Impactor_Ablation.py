"""
Script Name: Impactor_Ablation.py
Author: Emilia Vlahos
Date: 2024-10-20

Description:
    Simulates the dynamics and ablation of an impactor as it travels through a planetary atmosphere up to the (outermost) radiative-convective boundary.
    
    Inputs:
    - Atmospheric properties (Altitude, Temperature, Pressure, Density, Mean molecular weight, Adiabatic index).
        - Atmospheric profile is characterized by parameters of metallicity, opacity, disk density, orbital distance, core mass, advective depth, and GCR.
    - Impactor properties (Material, size)
    
    Outputs:
    - Impactor properties (Surface Temperature, Velocity, Radius, Mass/Mass Ablated) 
    - Corresponding atmospheric properties

"""

# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import pickle
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
# Import required Modules
from Constants import * 
from Functions import *

#%%
### Atmospheric Properties:
M_core = 5*M_e # planet core mass [kg]; MUST correspond with index d

## select desired atmospheric properties
# indices:
# m_index: metallicity; 0 to 6 for Z = 0.020, 0.002, 0.200, 0.400, 0.500, 0.800
# a: opacity; 0 to 1 for dusty and dust-free
# b: disk density; 0 to 2 for rho_disk = [1, 0.01, 0.0001]*rho_MMEN
# c: orbtial distance; 0 to 3 for a = [0.1, 0.3, 1.0, 3.5]*au
# d: core mass; 0 for M_core = 5 earth mass
# e: advective depth; 0 for alpha_adv = [0.3] 
# f: GCR value; 0 for GCR = 0.0591 (~6 %)

label, path = profile_label(0, 0, 1, 3, 0, 0, 0)
print(f"Running atmosphere with label:{label}")

#%% Load data from corresponding profile

input_dir = 'Users/emiliavlahos/atmospheric-pollution-repo/input_data' # input directory where necessary folders are stored
envelope_data = f"/{input_dir}/{path}/" # directory corresponding to specific profile

def load_and_convert(zone):
    """Loads data within a given atmospheric region and converts units to MKS
    Args:
    - Zone (str): dominant form of mixing describing the atmospheric region 
    
    Returns:
    - radius_values, T_values, rho_values, P_values (np.ndarray, floats): atmospheric data in MKS
    """
    
    directory = envelope_data + zone
    os.chdir(directory)
    
    radius_values = np.loadtxt('radius', delimiter=",")
    radius_values *= R_e # [earth radii] to [m]

    rho_values = np.loadtxt('density', delimiter=",")
    rho_values *=  1000 # [g/cm^3] to [kg/m^3]

    T_values = np.loadtxt('temp', delimiter=",")

    P_values = np.loadtxt('pressure', delimiter=",")
    P_values *= 0.1 # [dyne/cm^2] to [Pa]
    
    mu_values = np.loadtxt('molecular_weight', delimiter=",")
    molm_values = mu_values * 1.008e-3 # convert to molar masss in [kg/mol]
    
    gamma_values = np.loadtxt('ad_index', delimiter=",")
    
    return radius_values, rho_values, T_values, P_values, molm_values, gamma_values

# import data for full envelope (depending on the profile some zones may be empty yielding a user warning)
r_out, rho_out, T_out, P_out, molm_out, gamma_out = load_and_convert('outer/') # outer zone
r_rad, rho_rad, T_rad, P_rad, molm_rad, gamma_rad = load_and_convert('radiative/') # radiative zone
r_conv2, rho_conv2, T_conv2, P_conv2, molm_conv2, gamma_conv2 = load_and_convert('conv2/') # outer convective zone
r_conv1, rho_conv1, T_conv1, P_conv1, molm_conv1, gamma_conv1 = load_and_convert('conv1/') # inner convective zone

# outermost radius
initial_altitude = r_out[0]

# radiative-convective boundary; if the outer convective zone is non-zero, we take the outer rcb
if len(r_conv2)!=0: 
    rcb = r_conv2[0] 
else:
    rcb = r_conv1[0]

#%% Combine zones and arrange by distance from planet's center to ensure proper order
r = np.concatenate((r_out, r_rad, r_conv2, r_conv1))
rho = np.concatenate((rho_out, rho_rad, rho_conv2, rho_conv1))
T = np.concatenate((T_out, T_rad, T_conv2, T_conv1))
P = np.concatenate((P_out, P_rad, P_conv2, P_conv1))  
molm = np.concatenate((molm_out, molm_rad, molm_conv2, molm_conv1))
gamma= np.concatenate((gamma_out, gamma_rad, gamma_conv2, gamma_conv1))  

sorted_indices = np.argsort(r)
r = r[sorted_indices]
T = T[sorted_indices]
rho = rho[sorted_indices]
P = P[sorted_indices] 
molm = molm[sorted_indices] 
gamma = gamma[sorted_indices] 
    
#%% Interpolate temperature, pressure, density, molar mass, and adiabatic index as functions of radius
interp_T = interp1d(r, T, kind='linear')
interp_rho = interp1d(r, rho, kind='linear')
interp_P = interp1d(r, P, kind='linear')
interp_gamma = interp1d(r, gamma, kind='linear')
interp_molm = interp1d(r, molm, kind='linear')

#%% Impactor properties

# Select impactor species and its specific properties
material_label = 'Water' # Quartz, Water, or Carbon
SelectedMaterial = material_classes[material_label]
rho_i = SelectedMaterial(np.nan, np.nan).rho # impactor density of selected material [kg/m^3]

# Select general impactor properties
size_label = "10_m" 
initial_radius = 10 # initial impactor radius [m]

initial_mass =  rho_i * 4/3 * np.pi * initial_radius ** 3 # [kg]
initial_velocity = - np.sqrt(2 * G * M_core / initial_altitude) # (escape speed at initial altitude) [m/s]

#%% EOM of impactor entering atmosphere
initial_guess = 300 # initial guess for surface temperature of the impactor: I suggest 1000 @ 0.1 au and 300 @ 3.5 au
guess = initial_guess
# Function to define the system of ODEs
def ode(t, y):
    global guess  # Global variable to keep track of surface temperature guess across iterations
    
    # Unpack the current state variables from the input array y
    r = y[0]  # Atmospheric altitude [m]
    v = y[1]  # Impactor velocity [m/s]
    m = y[2]  # Impactor mass [kg]
    
    # Gas properties at the current altitude
    T = interp_T(r)      
    P = interp_P(r)       
    rho = interp_rho(r)   
    gamma = interp_gamma(r) 
    molm = interp_molm(r) 
    
    impactor = SelectedMaterial(v, m) # object representing the impactor with current velocity and mass
    gas = GasProperties(r, T, P, rho, gamma, molm) # object to represent the current atmospheric conditions
    
    # If Mach number is greater than 1, use shocked gas properties
    Ma = mach_number(impactor, gas)
    if Ma > 1:
        gas = ShockedGasProperties(r, T, P, rho, gamma, molm, Ma)
    
    T_s = surface_temperature(impactor, gas, guess)
    guess = T_s # Update the guess for the next iteration
    
    # Define the rates of change for the system (ODEs)
    drdt = v
    dvdt = - gravitational_acceleration(r, M_core) + drag_acceleration(T_s, impactor, gas)
    dmdt = - impactor.sublimation_rate(T_s, gas)
    
    return [drdt, dvdt, dmdt] # derivatives to be used in the ODE solver

#%% terminates solver at rad/conv boundary or if fully evapourated
def event1(t, y):
    # track mass to determine if it gets significantly small to determine that the impactor has fully evapourated
    return y[2] - 1e-7
def event2(t, y):
    # track radial distance to stop solving when the impactor is outside the current zone
    return y[0] - rcb

event1.terminal = True
event2.terminal = True

#%% Solve 
time_values = np.linspace(0, 10000000, 200000) # ajust as needed

# Use solve_ivp to numerically solve the equations of motion
# ajust error tolerance as needed; should be low for Carbon but can be raised for Quartz/Water
solution = solve_ivp(ode, [time_values[0], time_values[-1]], [initial_altitude, initial_velocity, initial_mass], t_eval=time_values, events=(event1, event2), rtol=1e-12, atol=1e-12)

t = solution.t
r = solution.y[0]  # Array of altitude values
v = solution.y[1]  # Array of velocity values
m = solution.y[2]  # Array of impactor mass values

# interpolate over the r to ensure an evenly spaced spacial grid
interp_v = interp1d(r, v, kind='linear')
interp_m = interp1d(r, m, kind='linear')

r = np.arange(r[0], r[-1], -100000)
v = interp_v(r)
m = interp_m(r)

# solve for the remaining parameters
T = interp_T(r)
P = interp_P(r)
rho = interp_rho(r)
molm = interp_molm(r)
gamma = interp_gamma(r)

#%%
# define impactor and gas objects to store data
impactor_data = SelectedMaterial(v, m)
gas_data = GasProperties(r, T, P, rho, gamma, molm)
gas_data = ShockedGasProperties(r, T, P, rho, gamma, molm, mach_number(impactor_data, gas_data))

#%% solve for impactor surface temeperature at each altitude index
# This is the 'bottle neck' step and can be skipped if desired once the surface temperature solver has been sufficiently tested at the relevant parameters
Ma = gas_data.Ma
guess = initial_guess
T_s = np.zeros_like(r) 
for i in range(len(r)):
    impactor = SelectedMaterial(v[i], m[i])
    gas = ShockedGasProperties(r[i], T[i], P[i], rho[i], gamma[i], molm[i], Ma[i])
    T_s[i] = surface_temperature(impactor, gas, guess) 
    guess = T_s[i] # update guess   
    
#%% Useful plot's to check that solver funtions correctly
fig = plt.figure()
plt.plot(r, radiative_power(T_s, impactor_data, gas_data), label='Radiative Power')
plt.plot(r, frictional_power(T_s, impactor_data, gas_data), label='Frictional Power')
plt.plot(r, evaporative_power(T_s, impactor_data, gas), label='Evaporative Power')
plt.plot(r, stepwise_temp(T_s, impactor_data, gas_data))
plt.xscale('log')
plt.yscale('symlog')
plt.xlabel('Altitude, m')
plt.ylabel('Power, W')
plt.legend()

fig = plt.figure()
plt.plot(r, gas_data.T, '-k', label='Gas Temperature')
plt.plot(r, gas_data.get_temperature(), '--k', label='Shocked Gas Temperature')
plt.plot(r, T_s, label='Impactor Surface Temperature')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Altitude, m')
plt.ylabel('Temperature, K')
plt.legend()

fig = plt.figure()
plt.plot(r, impactor_data.m, label='Impactor Mass')
plt.axvline(rcb)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Altitude, m')
plt.ylabel('Mass, kg')
plt.legend()

#%%% Save data 
interim_dir = '/Users/emiliavlahos/atmospheric-pollution-repo/interim_data/' # directory where intermidiate ablation data is stored
directory = f"{interim_dir}ablation_profiles/{material_label}/{size_label}/{path}"

os.makedirs(directory, exist_ok=True)
os.chdir(directory)

# Save gas data to a file
with open('gas_data.pkl', 'wb') as file:
    pickle.dump(gas_data, file)
    
# Save impactor data to a file
with open('impactor_data.pkl', 'wb') as file:
    pickle.dump(impactor_data, file)
    
# save impactor surface temperature data as a text file
np.savetxt('impactor_temperature', T_s, delimiter=",")
