## Atmospheric Pollution
This repository contains code I developed to explore the impact of solid accretion on atmospheric metallicity. The project focuses on modeling the ablation of solid impactors as 
they pass through an atmosphere and quantifying how much of the evaporated material is incorporated into the atmosphere.

The code is designed to be run for any general atmosphere, characterized by parameters of metallicity, opacity, disk density, orbital distance, core mass, advective depth, and 
gas-to-core ratio. For a given atmospheric profile, the code requires temperature, pressure, density, mean molecular weight, and adiabatic index data as a function of altitude, 
divided into advective, radiative, and convective zones.

To demonstrate its use, I’ve included example envelope data calculated by Vincent Savignac based on his model in Savignac & Lee (2024). 
If you use or adapt this code or any of the profiles provided, please cite [Vlahos et al (2024)](https://ui.adsabs.harvard.edu/abs/2024arXiv240913820V/abstract) 
and [Savignac & Lee (2024)](https://ui.adsabs.harvard.edu/abs/2024ApJ...973...85S/abstract).

### Description of Each Script

- **Constants.py**: Defines constants and set values used throughout the project, which are imported into the other scripts.
  
- **Functions.py**: A collection of classes and functions that are either used across scripts or provide intermediate steps in the calculations.
  
- **Impactor_Ablation.py**: Calculates the ablation of a single solid impactor within a given atmosphere. The material of the impactor can be selected
via the `material_label` to model the ablation of quartz, water ice, or refractory carbons. Atmospheric parameters are chosen independently, but they should
correspond to the expected solid disk composition (e.g., water ice should only be selected for atmospheres outside the snow line). The initial guess for the surface
temperature solver should also reflect the ambient temperature, which changes with orbital distance.
  
- **Mass_Fraction.py**: Calculates and plots the post-accretion atmospheric mass fraction from the ablation data for a specified atmosphere and impactor.
  Accretion can be modeled either in situ or migrating, and both solid disk density and accretion timescale can be adjusted to modify the net amount of accretion.
  
- **Carbon_to_Oxygen.py**: Calculates and plots the post-accretion atmospheric carbon-to-oxygen ratio from the ablation data of water ice and refractory carbon into
  a specified atmosphere. Accretion can similarly be modeled either in situ or migrating, with the solid disk density and accretion timescale adjustable to modify the net amount of accretion.
