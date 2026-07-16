# bomb-calorimeter-analysis
Python-based analysis of bomb calorimeter experimental data.
# Bomb Calorimeter Data Analysis

This repository contains a Python-based analysis tool for bomb calorimeter experimental data. It processes raw temperature-time data from an Isoperibol Calorimeter (simulated IC210) to calculate the heat of combustion of hydrocarbon fuels, adhering to standard laboratory practices.

## Overview
In experimental physics and thermodynamics, accurately determining the heat of combustion involves meticulous data processing. This project automates the workflow:
1. Parsing raw temporal and thermal data.
2. Visualizing the temperature change during combustion.
3. Calculating the heat of combustion ($H$) using the relation $Q = C_{cal}\Delta T$ and $H = \frac{Q}{m}$.
4. Performing uncertainty and error propagation:
$$
\Delta Q = Q\sqrt{\left(\frac{\Delta C}{C}\right)^2 + \left(\frac{\Delta T}{\Delta T}\right)^2}
$$

## Supported Standards
The calculation modules are designed with the infrastructure to comply with major international calorimetry standards, including:
* **Liquid Fuels:** ASTM D240, ASTM D4809
* **Solid Mineral Fuels & Coal:** ISO 1928, ASTM D5865
* **Biomass & Biofuels:** ISO 18125
* **Waste & Refuse-Derived Fuels:** ASTM D5468, ASTM E711

## Tech Stack
* Python 3
* Pandas (Data I/O and manipulation)
* Matplotlib (Visualization)
* SciPy (Optimization and curve fitting)

## Results & Error Analysis
The thermodynamic calculations incorporate strict error propagation based on GUM (Guide to the Expression of Uncertainty in Measurement) standards. 

By accounting for the instrumental uncertainties of the temperature sensor (± 0.0005 °C), the analytical balance (± 0.0001 g), and the inherent calibration variance, the specific heat of combustion is determined with high precision:

* **Final Result:** H = 26477.7 ± 36.8 J/g
* **Relative Error:** 0.139 %

Achieving a relative uncertainty of less than 0.2% demonstrates both the reliability of the simulated IC210 system data and the robustness of the automated processing pipeline.
