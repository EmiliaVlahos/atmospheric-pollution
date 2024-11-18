"""
Script Name: Constants.py
Author: Emilia Vlahos
Date: 2024-10-20
Description: Stores constants used accross scripts
"""

# Astronomical contants
G = 6.6743e-11 # Universal gravitational constant [N*m^2/kg^2]
au = 1.496e+11 # Astronomical unit [m]
M_e = 5.9722e24 # Earth's mass [kg]
R_e = 6.3781e6 # Earth's radius [m]
M_s = 1.989e30 # Solar mass [kg]

# Stat-mech constants
sigma = 5.67037e-8 # Stefan-Boltzmann constant [W/m^2*K^4]
R = 8.3145 # Universal gas constant [J/mol*K]
k_b = 1.380649e-23 # Boltzmann constant [kg*m^2/s^2*K]

# Properties of Hydrogen and Helium
molm_H2 = 2.016e-3 # Molar mass of H2 [kg/mol]
molm_He = 4.003e-3 # Molar mass of He [kg/mol]
mm_H2 = 3.34e-27 # Molecular mass of H2 [kg]
mm_He = 6.65e-27 # Molecular mass of He [kg]
d_H2 = 2.89e-10 # Diameter of H2 molecule [m]
d_He = 2.60e-10 # Diameter of He atom [m]

# Properties of SiO2
molm_SiO2 = 6.008e-2 # Molar mass of SiO2 [kg/mol]
rho_SiO2 = 2.65e3 # Density of SiO2 [kg/m^3]
E_vap_SiO2 = 8.08e6 # Latent heat of vaporization of SiO2 [J/kg]

# Properties of H2O
molm_H2O = 1.8015e-2 # Molar mass of H2O [kg/mol]
rho_H2O = 917 # Density of H2O (ice) [kg/m^3]
E_vap_H2O = 2.83e6 # Latent heat of vaporization of H2O [J/kg]

# Properties of refractory carbons
rho_refcarb = 1.07e3 # Density of type I kerogen [kg/m^3]
E_vap_refcarb = 3.2e6 # Latent heat of sublimation of carbonous chondrites [J/kg]
# based on assummed reftractory carbon composition of C100H77O14N3.5S2
C_per_kg = 61.96 # moles of Carbon per kg of refractory carbon
O_per_kg = 8.674 # moles of Oxygen per kg of refractory carbon

# Properties of disk 
molm_disk = 2.37 * 1.008e-3 # Molar mass of disk gas [kg/m^3]

# Solar abundances from Grevesse & Noels (1993)
X = 0.7 # Solar Hydrogen mass fraction
Y = 0.28 # Solar Helium mass fraction
Z = 0.02 # Solar metallic species mass fraction
## metallic species mass fraction 
# with respect to net metallicity Z
silicon_mf = 0.040536 
oxygen_mf = 0.482398 
carbon_mf = 0.173312
# with respect to total mass
solar_silicon = silicon_mf * Z 
solar_oxygen = oxygen_mf * Z 
solar_carbon = carbon_mf * Z
# Carbon and Oxygen
oxygen_nf = 0.513001 # Oxygen number fraction with respect to net metallicity Z
carbon_nf = 0.245538 # Carbon number fraction with respect to net metallicity Z
solar_CtoO = carbon_nf / oxygen_nf # solar carbon to oxygen number ratio
solid_CtoO = solar_CtoO / 2 # carbon to oxygen ratio within the solid disk
molm_O = 1.5999e-2 # Molar mass of oxygen atom [kg/mol]
molm_C = 1.2011e-2 # Molar mass of carbon atom [kg/mol]
water_disk_fraction = 0.823 # Mass fraction of water ice within solid disk
carbon_disk_fraction = 0.177 # Mass fraction of refractory carbons within solid disk