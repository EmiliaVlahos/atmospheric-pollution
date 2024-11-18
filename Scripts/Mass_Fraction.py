"""
Script Name: Mass_Fraction.py
Author: Emilia Vlahos
Date: 2024-10-27

Description:
    Calculates and plots the post-accretion atmospheric mass fraction of accreted materials 
    
    Inputs:
    - Atmospheric (planet) properties (Metallicity, Opacity, Disk density, Orbital distance, Core mass, Advective depth, GCR).
    - Impactor properties (Material, radius)
    - Solid disk density (fraction of MMEN)
    - Formation type (in situ, migrating)
    - Accretion timescale

"""

# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os

# Import required Modules
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from Constants import * 
from Functions import *

# directory from which data is imported
interim_dir = '/Users/emiliavlahos/Atmospheric_Pollution_repo/interim_data/' # directory where intermidiate ablation data is stored

#%% Select atmospheric properties
# indices:
# m_index: metallicity; 0 to 6 for Z = 0.020, 0.002, 0.200, 0.400, 0.500, 0.800
# a: opacity; 0 to 1 for dusty and dust-free
# b: disk density; 0 to 2 for rho_disk = [1, 0.01, 0.0001]*rho_MMEN
# c: orbtial distance; 0 to 3 for a = [0.1, 0.3, 1.0, 3.5]*au
# d: core mass; 0 for M_core = 5 earth mass
# e: advective depth; 0 for alpha_adv = [0.3] 
# f: GCR value; 0 for GCR = 0.0591 (~6 %)
label, path = profile_label(0, 0, 1, 3, 0, 0, 0)

f_gas =  0.01 # Gaseous disk density as a fraction of MMEN; MUST correspond with index b
orbital_distance = 3.5*au # (Initial) orbital distance [m]; MUST correspond with index c
M_core = 5*M_e # planet core mass [kg]; MUST correspond with index d

#%% select impactor properties
species_label = 'Water' # impactor species; Quartz, Water, or Carbon
size_label = "10_m" # initial impactor radius; eg. 10_m

directory = f"{interim_dir}ablation_profiles/{species_label}/{size_label}/{path}"
os.chdir(directory)

# Load gas data
with open('gas_data.pkl', 'rb') as file:
    gas_data = pickle.load(file)
# Load impactor data
with open('impactor_data.pkl', 'rb') as file:
    impactor_data = pickle.load(file)

#%% select other properties
solid_disk_fractions = [1e-6, 1e-5, 1e-4] # solid disk densities as a fraction of MMEN
timescale = 2e6*3.154e+7 # accretion timescale [seconds]
formation_type = 'in situ' # 'in situ' or 'migrating'

#%% plot mass fraction for a range of solid disk densities
colors = sns.color_palette("cubehelix", 4).as_hex()[::-1]

# set up plot
fig = plt.figure()
plt.gca().invert_xaxis() 
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Gas Pressure (bar)')
plt.ylabel(f'{species_label} Mass Fraction')

# calulate and plot mass fraction for a range of solid disk densities
for i, f_solid in enumerate(solid_disk_fractions):
    mf = mass_fraction(impactor_data, gas_data, orbital_distance, f_solid, f_gas, M_core, timescale, formation_type) # calculate mass fraction
    
    # plot mass fraction relative to gas pressure (converted to bar)
    plt.plot(gas_data.P * 1e-5, mf, color=colors[i], label=f"$f_{{\mathrm{{solid}}}} = {f_solid}$") 
    
plt.plot(gas_data.P * 1e-5, max_fraction(impactor_data, gas_data), color=colors[3], label='Saturation Limit') 
plt.axhline(impactor_data.solar, color='k', ls='dotted')
plt.legend()
