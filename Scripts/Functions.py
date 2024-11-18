"""
Script Name: Functions.py
Author: Emilia Vlahos
Date: 2024-10-20
Description: Functions and classes used accross scripts
"""

import numpy as np
from scipy.optimize import root
from Constants import *

#%% Function to navigate atmosphere parameter space

# Select envelope parameters and generate corresponding label to retrieve data
metallicity_names = ['Z=0.020','Z=0.002','Z=0.200','Z=0.400','Z=0.500','Z=0.800'] # metallicity (value for Z)
opacities = ['Dusty','Dust-free'] # opacity scheme
disk_names = ['1_MMEN','0.01_MMEN','0.0001_MMEN'] # disk density depletion from MMEN
distance_names=['0.1_AU','0.3_AU','1_AU','3.5_AU'] # orbital distance
core_names = ['5_Mearth'] # core mass
adv_names = ['0.3_adv'] # advection depth
gcr_names = ['GCR=0.0591'] # gas to core ratio

def profile_label(m, a, b, c, d, e, f):
    """ Generate a profile label based on the given indices.
    Args: Indices (int)
    - m: metallicity; 0 to 6 for Z = 0.020, 0.002, 0.200, 0.400, 0.500, 0.800
    - a: opacity; 0 to 1 for dusty and dust-free
    - b: disk density; 0 to 2 for rho_disk = [1, 0.01, 0.0001]*rho_MMEN
    - c: orbtial distance; 0 to 3 for a = [0.1, 0.3, 1.0, 3.5] au
    - d: core mass; 0 for M_core = 5 earth mass
    - e: advective depth; 0 for alpha_adv = [0.3] 
    - f: GCR value; 0 for GCR = 0.0591 (~6 %)

    Returns:
    - label (str): A formatted label string.
    - path (str): A directory path to retrieve the specific profile from the data folder
    """
    
    label = f"{metallicity_names[m]}_{opacities[a]}_{disk_names[b]}_{distance_names[c]}_{core_names[d]}_{adv_names[e]}_{gcr_names[f]}"
    path = f"{metallicity_names[m]}/{opacities[a]}/{disk_names[b]}/{distance_names[c]}/{core_names[d]}/{adv_names[e]}/{gcr_names[f]}"
    
    return label, path

#%% Classes to store gas properties
class GasProperties:
    """
    Class representing the thermophysical properties of a gas.

    Attributes:
        r (float or array): Altitude [m].
        T (float): Temperature [K].
        P (float): Pressure [Pa].
        rho (float): Density [kg/m^3].
        gamma (float): Adiabatic index (dimensionless).
        molm (float): Molar mass [kg/mol].
    """
    
    def __init__(self, r, T, P, rho, gamma, molm):
        self.r = r # altitude [m]
        self.T = T # temperature [K]
        self.P = P # pressure [Pa]
        self.rho = rho # density [kg/m^3]
        self.gamma = gamma # adiabatic index (dimensionless)
        self.molm = molm # molar mass [kg/mol]
        
    def get_temperature(self):
        return self.T

    def get_pressure(self):
        return self.P

    def get_density(self):
        return self.rho
    
    def specific_heat(self):
        """Returns specific heat of the gas [J/kg*K]."""
        T = self.get_temperature()
        P = self.get_pressure()
        rho = self.get_density()
        grad_ad = (self.gamma - 1) / self.gamma
        return P / (T * rho * grad_ad) # [J/kg*K]

    def dyn_vis(self):
        """Returns dynamic viscosity of hydrogen at gas conditions [Pa*s]."""
        T = self.get_temperature()
        return (5 / (16 * d_H2 ** 2)) * np.sqrt(mm_H2 * k_b * T / np.pi) # [Pa*s]

    def kin_vis(self):
        """Returns kinematic viscosity of hydrogen at gas conditions [m^2/s]."""
        rho = self.get_density()
        eta = self.dyn_vis()
        return eta / rho # [m^2/s]

    def therm_cond(self):
        """Returns thermal conductivity of hydrogen at gas temperature [W/m*K]."""
        T = self.get_temperature()
        return (1434 + 1.257 * (T - 1200)) * 1e-6 * 418.4 # [W/m*K]

    def sound_speed(self):
        """Returns sound speed in the gas [m/s]."""
        T = self.get_temperature()
        return np.sqrt(self.gamma * R * T / self.molm) # [m/s]
    
    def envelope_mass(self):
        """"Calculates the mass of the atmospheric gas contained in layer of envelope [kg]. 
        - not applicable for local instances of GasProperties"""
        # Check if self.r is an array.
        if not isinstance(self.r, (np.ndarray, list)):
            raise ValueError("The altitude (r) must be an array or a list.")

        shell_volume = np.append(((4/3) * np.pi * abs(np.diff(self.r**3))), np.nan) # append nan value to keep consistent length
        return self.rho * shell_volume # [kg]
        
        
class ShockedGasProperties(GasProperties):
    """
    Sub-class of gas properties; calculates and utalizes shocked-gas properties.
    """
    
    def __init__(self, r, T, P, rho, gamma, molm, Ma):
        super().__init__(r, T, P, rho, gamma, molm) # includes pre-shock temperature, pressure, and density
        self.Ma = Ma # pre-shock mach number
        
        self.T_shocked = self._calculate_shocked_temperature()
        self.P_shocked = self._calculate_shocked_pressure()
        self.rho_shocked = self._calculate_shocked_density()
        
    def _calculate_shocked_temperature(self):
        """Returns the shocked gas temperature [K]."""
        Ma = np.maximum(self.Ma, 1) # if Mach number is less than one, return pre-shock temperature
        factor = (2 * self.gamma * Ma**2 - (self.gamma-1)) * (2 + (self.gamma-1) * Ma**2) / ((self.gamma+1)**2 * Ma**2)
        return self.T * factor # [K]
        
    def _calculate_shocked_pressure(self):
        """Returns the shocked gas pressure [Pa]."""
        Ma = np.maximum(self.Ma, 1) # if Mach number is less than one, return pre-shock pressure
        factor = (2 * self.gamma * Ma**2 - (self.gamma-1))/ (self.gamma+1)
        return self.P * factor # [Pa]

    def _calculate_shocked_density(self):
        """Returns the shocked gas density [kg/m^3]."""
        Ma = np.maximum(self.Ma, 1) # if Mach number is less than one, return pre-shock density
        factor = (self.gamma+1) * Ma**2 / (2 + (self.gamma-1)*Ma**2)
        return self.rho * factor # [kg/m^3]
    
    def get_temperature(self):
        return self.T_shocked

    def get_pressure(self):
        return self.P_shocked

    def get_density(self):
        return self.rho_shocked
    
    
#%% Classes to store impactor properties  
class ImpactorProperties:
    """
    Class representing the physical and material properties of an impactor.

    Attributes:
        v (float): Impactor velocity in meters per second [m/s].
        m (float): Impactor mass in kilograms [kg].
    """
    
    def __init__(self, v, m):
        self.v = v # impactor velocity [m/s]
        self.m = m # impactor mass [kg]
        
    alpha = 1 # impactor absortivity
          
    def radius(self):
        """Returns the impactor's radius as a function of it's mass and density [m]."""
        return np.cbrt(3 * self.m / (4 * np.pi * self.rho)) # [m]
        
    def sublimation_rate(self, T_s):
        """Returns the impactor's mass-loss rate specific to its material [kg/s]."""
        raise NotImplementedError("This method should be overridden in a subclass.")
        
    def mass_lost(self):
        """"Calculates the mass lost by the impactor in each atmospheric layer [kg]. 
        - Not applicable for instantaneous instances of ImpactorProperties (mass must evolve)."""
        # Check if self.m is an array.
        if not isinstance(self.m, (np.ndarray, list)):
            raise ValueError("The mass (m) must be an array or a list.")
        
        return np.append(abs(np.diff(self.m)), 0) # append zero to keep consistent length
        
    
class QuartzProperties(ImpactorProperties):
    """
    Subclass representing SiO2 impactors and their ablation processes.
    """
    def __init__(self, v, m):
        super().__init__(v, m)
        
    molm = molm_SiO2 # [kg/mol]
    rho = rho_SiO2 # [kg/m^3]
    E_vap = E_vap_SiO2 # [J/kg]
    solar = solar_silicon # solar silicon mass fraction for comparison
    
    def vapor_pressure(self, T_s):
        """Returns the temperature dependant vapor pressure of SiO2 [Pa]."""
        return  0.1 * np.exp(29.5 - 46071.4/(T_s + 58.9)) #[Pa]
    
    def sublimation_rate(self, T_s, gas=None):
        """Returns the (absolute) mass loss rate of an SiO2 impator as a function of its surface temperature [kg/s]."""
        r_i = self.radius()
        P_vap = self.vapor_pressure(T_s)
        return 4 * np.pi * r_i**2 * P_vap * np.sqrt(self.molm / (2 * np.pi * R * T_s)) # [kg/s]
    
    def max_deposition(self, gas):
        """Returns the maximum mass that the impactor can deposite into an atmosphere (with properties gas) before it condenses out."""
        P_vap = self.vapor_pressure(gas.T) # assuming evaporated material is heated to the local gas temperature
        f_saturation = self.molm * P_vap / (gas.molm * gas.P)
        return gas.envelope_mass() * f_saturation
    

class WaterProperties(ImpactorProperties):
    """
    Subclass representing H2O impactors and their ablation processes.
    """
    molm = molm_H2O # [kg/mol]
    rho = rho_H2O # [kg/m^3]
    E_vap = E_vap_H2O # [J/kg]
    solar = solar_oxygen # solar oxygen mass fraction for comparison
    
    def vapor_pressure(self, T_s):
        """Returns the temperature dependant vapor pressure of H2O [Pa]."""
        
        T_m = 272.84 # Melting point of H2O [K]
        T_cr = 647.096 # Liquid-vapor critical temperature of H2O [K]
        P_cr = 2.2064e8 # Vapor pressure at critical temperature [dyne/cm^2]
        
        def low_Ts(T_s):
            # Returns the vapor pressure of H2O below 272 K.
            a0 = -2445.5646
            a1 = -3.63
            a2 = -0.01677006
            a3 = 1.20514e-5
            a4 = 8.2312
            aux = (a0 / T_s) + a1 + (a2 * T_s) + (a3 * T_s**2)  + (a4 * np.log10(T_s)) # log of vapor pressure in [dyne/cm^2]
            return 0.1 * np.power(10, aux) # reverse log and convert [dyne/cm^2] to [Pa]

        def mid_Ts(T_s):
            # Returns the vapor pressure  of H2O between 272 K and 647 K.
            a0 = -7.85951783
            a1 = 1.84408259
            a2 = -11.7866497
            a3 = 22.6807411
            a4 = -15.9618719
            a5 = 1.80122502
            theta_s = 1 - T_s / T_cr
            aux = (T_cr / T_s) * (a0 * theta_s + a1 * theta_s**(3/2) + a2 * theta_s**3 + 
                                  a3 * theta_s**(7/2) + a4 * theta_s**4 + a5 * theta_s**(15/2))  
            return 0.1 * P_cr * np.exp(aux) # reverse log and convert [dyne/cm^2] to [Pa]
        
        # Apply peicewise logic
        P_vap = np.zeros_like(T_s, dtype=float)
        P_vap[T_s <= T_m] = low_Ts(T_s[T_s <= T_m])
        P_vap[T_s > T_m] = mid_Ts(T_s[T_s > T_m])
        P_vap[T_s > T_cr] = np.nan # vapor pressure model fails if impactor surface temp surpasses T_cr (can't have solid ice in very hot regions)
        return P_vap # [Pa]
    
    def sublimation_rate(self, T_s, gas=None):
        """Returns the (absolute) mass loss rate of an H2O impator [kg/s]."""
        r_i = self.radius()
        P_vap = self.vapor_pressure(T_s)
        return 4 * np.pi * r_i**2 * P_vap * np.sqrt(self.molm / (2 * np.pi * R * T_s)) # [kg/s]
    
    def max_deposition(self, gas):
        """Returns the maximum mass that the impactor can deposite into an atmosphere (with properties of gas) before it condenses out [kg]."""
        P_vap = self.vapor_pressure(gas.T) # assuming evaporated material is heated to the local gas temperature
        f_saturation = self.molm * P_vap / (gas.molm * gas.P)
        intermidiate = gas.envelope_mass() * f_saturation
        return  np.where(np.isnan(intermidiate), np.inf, intermidiate) # infinity when vapor pressure isn't defined since deposition is not unlimited here
    
    
class CarbonProperties(ImpactorProperties):
    """
    Subclass representing carbonous impactors and their ablation processes.
    """
    rho = rho_refcarb # [kg/m^3]
    E_vap = E_vap_refcarb # [J/kg]
    solar = np.nan # solar carbon is not a great means of comparison since composition of refractory carbons is more complex
    
    def impactor_conductivity(self, T_s):
        """Returns the temperature dependant thermal conductivity [W/mK] of kerogen rich oil shales."""
        T_C = T_s - 273.15 # convert Kelvin to Celsius
        return 1.139 - 6.29e-4*T_C # [W/mK]
    
    def isothermal_mass(self, T_s, gas):
        """Returns the total mass [kg] of an isothermal outer layer heated to T_s."""
        T = gas.get_temperature()
        k_carb = self.impactor_conductivity(T_s)
        r_i = self.radius()
        d_iso = np.minimum((0.3 * k_carb / (sigma * abs(T_s - T)**3)), r_i) # depth of isothermal layer
        return (4/3) * np.pi * self.rho * (r_i**3 - (r_i - d_iso)**3) # [kg]
        
    def sublimation_rate(self, T_s, gas):
        """Returns the (absolute) mass loss rate of a carbonous impator [kg/s]."""
        m_iso = self.isothermal_mass(T_s, gas) 
        # Kinetic rate parameters of Kerogen type I
        E_a = 2.31e05 # activation energy [J/mol]
        A = 1.70e14 # frequency cst [s^-1]
        k = A * np.exp(-E_a / (R*T_s)) # rate constant [s^-1]
        return m_iso * k # [kg/s]
        
    def max_deposition(self, gas):
        """Returns an array of np.inf to keep code consistent between materials."""
        return np.full_like(gas.r, np.inf)

# Map the material type to the corresponding class
material_classes = {
    'Quartz': QuartzProperties,
    'Water': WaterProperties,
    'Carbon': CarbonProperties
}


#%% function of the impacor / gas instances utalized in ablation calculations

def reynolds_number(impactor, gas):
    """Returns Reynolds number of impactor's flow through the gas."""
    nu = gas.kin_vis()
    r_i = impactor.radius()
    return abs(impactor.v) * 2 * r_i / nu

def mach_number(impactor, gas):
    """Returns Mach number of impactor's flow through the gas."""
    cs = gas.sound_speed()
    return abs(impactor.v) / cs

def drag_coefficient(T_s, impactor, gas):
    """Returns the drag coefficient of the impactor's flow through the gas."""
    T = gas.get_temperature()
    Re = reynolds_number(impactor, gas)
    Ma = mach_number(impactor, gas)
    Kn = Ma / Re   
    
    C_1 = (24 / Re) * (1 + 0.15 * Re ** 0.678) + 0.407 * Re / (Re + 8710)
    C_2 =  10 ** ((2.5 * (Re / 312) ** 0.6688)/(1 + (Re / 312) ** 0.6688))
    H = 4.6 /(1 + Ma) + 1.7 * np.sqrt(T_s / T)
    C_d = 2 + (C_1 - 2) * np.exp(-3.07 * np.sqrt(gas.gamma) * Kn * C_2) + H * (1 / (np.sqrt(gas.gamma) * Ma)) * np.exp(-1 / (2 * Kn))
    C_d = np.array(C_d)
    C_d[Re > 3e5] = 0.2 # limit drag coefficient in turbulent regime
    return C_d

def radiative_power(T_s, impactor, gas):
    """Returns the power [W] radiated into the impactor from the gas."""
    T = gas.get_temperature()
    r_i = impactor.radius()
    return 4 * np.pi * r_i**2 * impactor.alpha * sigma * (T**4 - T_s**4) # [W]

def frictional_power(T_s, impactor, gas): 
    """Returns the power [W] associated with the impactors heating due to drag forces."""
    rho = gas.get_density()
    r_i = impactor.radius()
    cp_g = gas.specific_heat()
    k_g = gas.therm_cond()
    eta = gas.dyn_vis()
    Re = reynolds_number(impactor, gas)
    Ma = mach_number(impactor, gas)
    C_d = drag_coefficient(T_s, impactor, gas)

    Pr = eta * cp_g / k_g
    r_ = Pr ** (1/3)
    Nu_c = 2 + 0.459 * (Re)**0.55 * r_
    M_ = Ma / (0.428 * Ma * (gas.gamma + 1) / gas.gamma)
    Nu = Nu_c / (1 + 3.42 * M_ * Nu_c / (Re * Pr))
    f = (8 / gas.gamma) * (Nu / (Re * Pr)) * (r_ / C_d) # factor of kinetic energy converted to heat
    F_drag = (1/2) * np.pi * r_i**2 * rho * impactor.v**2 * C_d # drag force
    return f * F_drag * abs(impactor.v) # [W]

def evaporative_power(T_s, impactor, gas):
    """Returns the power [W] associated with the impactors cooling due to vaporization."""
    dm = - impactor.sublimation_rate(T_s, gas)
    return dm * impactor.E_vap # [W]

def stepwise_temp(T_s, impactor, gas):
    """Returns the sum of all power terms. A root finder will extract the impactor surface temperature by equating to zero."""
    P_rad = radiative_power(T_s, impactor, gas)
    P_fric = frictional_power(T_s, impactor, gas)
    P_evap = evaporative_power(T_s, impactor, gas)
    return P_rad + P_fric + P_evap # should be zero

def surface_temperature(impactor, gas, guess):
    """Returns the impactors surfaces temperature [K] by ballancing power into and out of the impactor."""
    sol = root(stepwise_temp, guess, args=(impactor, gas), method='hybr', tol=1e-12) # ajust tolerance as needed
    return sol.x[0]

def gravitational_acceleration(r, M_core):
    """Returns gravitational acceleration [m/s^2] as a function of r: distance to planet's center"""
    return G * M_core / (r**2) # [m/s^2]

def drag_acceleration(T_s, impactor, gas):
    """Returns impactor acceleration [m/s^2] due to drag."""
    rho = gas.get_density()
    r_i = impactor.radius()
    C_d = drag_coefficient(T_s, impactor, gas)
    return 1/2 * np.pi * r_i**2 * rho * impactor.v**2 * C_d / impactor.m # [m/s^2]

#%% function of the impacor / gas instances utalized in mass fraction calculations

def solid_surface_density(a, f_solid):
    """Returns the surface density [kg/m^2] of solid materials within the disk at an orbital distance a."""
    return f_solid * 6.2e3 * (a / (0.2 * au))**(-1.6) # [kg/m^2]

def orbital_frequency(a):
    """Returns keplerian frequency in [s^-1] at a given orbital distance a."""
    return np.sqrt(G * M_s / a**3) # [s^-1]

def accretion_rate(R_acc, a, f_solid):
    """Returns the solid accretion rate in [kg/s] into the planet from its outer radius R_acc."""
    Sigma_s = solid_surface_density(a, f_solid)
    Omega = orbital_frequency(a)
    return np.sqrt(np.pi/2) * R_acc**2 * Omega * Sigma_s # [kg/s]

def max_fraction(impactor, gas):
    """Returns the maximum (saturatecd) mass fraction of impactor mass relative to the total atmospheric mass (impactor + initial gas)"""
    m_max = impactor.max_deposition(gas) # max mass that can be held in each atm layer (if np.inf then unlimited) [kg]
    m_gas = gas.envelope_mass() # get corresponding atmospheric mass
    return m_max / (m_max + m_gas)

def total_solid_mass_insitu(impactor, gas, a, f_solid, timescale):
    """Returns the mass profile [kg] of the net impactor material held in the atmosphere after accretion over a time of timsecale for a planet accreting in situ."""
    R_acc = gas.r[0] # altitude at which we begin accretion calculations
    m_initial = impactor.m[0] # impactor mass before ablation
    dIdt = accretion_rate(R_acc, a, f_solid) / m_initial # [impactors/s]
    dm = impactor.mass_lost() 
    
    total_accreted_mass = dm * dIdt * timescale # total mass accreted into each atm layer [kg]
    m_max = impactor.max_deposition(gas) # max mass that can be held in each atm layer (if np.inf then unlimited) [kg]
    return np.minimum(total_accreted_mass, m_max) # [kg]

def OrbitalDistance(t, a_0, f_gas, M_core):
    """ Returns the time dependant orbital distance [m] of a migrating planet inwards from a_0 (derived by analytically integrating the planet's EOM)"""
    cst = 1.347224e3 * au**(41/35) * G**2 * M_core * f_gas * molm_disk/ (R * np.sqrt(G * M_s))
    return (a_0**(117/70) - (117/70)*cst*t)**(70/117) # [m]

def total_solid_mass_migrating(impactor, gas, a_0, f_solid, f_gas, M_core, timescale):
    """Returns the mass profile [kg] of the net impactor material held in the atmosphere after accretion over a time of timsecale for a planet accreting while migrating."""
    R_acc = gas.r[0] # altitude at which we begin accretion calculations
    m_initial = impactor.m[0] # impactor mass before ablation 
    dm = impactor.mass_lost() 
    
    time_values, dt = np.linspace(0, timescale, 200000, retstep=True)  # divide formation timescale into discrete steps
    total_accreted_mass = np.zeros_like(dm) # initialize array
    for time in time_values:
        
        a = OrbitalDistance(time, a_0, f_gas, M_core)  # update the planet's position
        
        # check that the planet is still outside the iceline
        if a < 1 * au:
            # if not, break loop
            t_iceline = time
            print('Ice accretion stops at:', t_iceline/3.154e+7, 'years')
            break
        
        dIdt = accretion_rate(R_acc, a, f_solid) / m_initial # accretion rate at a [impactors/s]
        total_accreted_mass += (dm * dIdt * dt) # total mass accreted into each atm layer over the timestep [kg]
        
    m_max = impactor.max_deposition(gas) # max mass that can be held in each atm layer (if np.inf then unlimited) [kg]
    return np.minimum(total_accreted_mass, m_max) # [kg]

def mass_fraction(impactor, gas, a, f_solid, f_gas, M_core, timescale, formation_type):
    """Returns the mass fraction of the net impactor mass relative to the total atmospheric mass (impactor + initial gas)"""

    # calculate mass for desired formation type: in situ or migrating
    if formation_type == 'in situ':
        m_solid = total_solid_mass_insitu(impactor, gas, a, f_solid, timescale)
    elif formation_type == 'migrating':
        m_solid = total_solid_mass_migrating(impactor, gas, a, f_solid, f_gas, M_core, timescale)
    else:
        raise ValueError("Invalid formation_type. Expected 'in situ' or 'migrating', got '{}'.".format(formation_type))
        
    m_gas = gas.envelope_mass() # get corresponding atmospheric mass
    
    return m_solid / (m_solid + m_gas)

#%%

def align_sizes(arr1, arr2, arr3):
    """ Trims the ends of three imput arrays so that they are compatible sizes."""
    min_length = min(len(arr1), len(arr2), len(arr3)) # Find the length of the shortest array
    
    # Trim both arrays to the length of the shorter one
    arr1_trimmed = arr1[:min_length]
    arr2_trimmed = arr2[:min_length]
    arr3_trimmed = arr3[:min_length]
    
    return arr1_trimmed, arr2_trimmed, arr3_trimmed

def carbon_to_oxygen(water_impactor, carbon_impactor, water_gas, carbon_gas, a, f_solid, f_gas, M_core, timescale, formation_type):
    """Returns the C/O ratio for a given atmosphere from the corresponding Water and Gas impactor/gas data objects."""

    # calculate mass for desired formation type: in situ or migrating
    if formation_type == 'in situ':
        m_water = water_disk_fraction * total_solid_mass_insitu(water_impactor, water_gas, a, f_solid, timescale)
        m_carbon = carbon_disk_fraction * total_solid_mass_insitu(carbon_impactor, carbon_gas, a, f_solid, timescale)
    elif formation_type == 'migrating':
        m_water = water_disk_fraction * total_solid_mass_migrating(water_impactor, water_gas, a, f_solid, f_gas, M_core, timescale)
        m_carbon = carbon_disk_fraction * total_solid_mass_migrating(carbon_impactor, carbon_gas, a, f_solid, f_gas, M_core, timescale)
    else:
        raise ValueError("Invalid formation_type. Expected 'in situ' or 'migrating', got '{}'.".format(formation_type))
        
    # Trim the mass deposition arrays and one of the corresponding atm gas mass arrays to ensure they are the same length.
    # This preserves alignment of atmospheric positions, as the grid is partitioned based on radius.
    m_water_trimmed, m_carbon_trimmed, m_gas = align_sizes(m_water, m_carbon, water_gas.envelope_mass())
    
    mol_oxygen = solar_oxygen * m_gas / molm_O # initial oxygen abundance (moles) assuming solar with respect to hydrogen
    mol_carbon = solar_carbon * m_gas / molm_C # initial carbon abundance (moles) assuming solar with respect to hydrogen
     
    mol_oxygen += m_water_trimmed/molm_H2O # additional  moles of oxygen deposited by icy impactors
    mol_oxygen += m_carbon_trimmed * O_per_kg # additional  moles of oxygen deposited by refractory carbons
    mol_carbon += m_carbon_trimmed * C_per_kg # additional  moles of carbon deposited by refractory carbons
    
    return mol_carbon / mol_oxygen
    