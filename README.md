# bomb-calorimeter-analysis

Python-based analysis of bomb calorimeter experimental data.

## Overview

This repository contains a Python-based analysis tool for bomb calorimeter experimental data. It processes raw temperature-time data from an Isoperibol Calorimeter (simulated IC210) to calculate the heat of combustion of hydrocarbon fuels, following standard laboratory practices.

In experimental physics and thermodynamics, accurately determining the heat of combustion involves meticulous data processing. This project automates the workflow:

1. Parsing raw temporal and thermal data (including instrument metadata such as sample mass and test standard).
2. Visualizing the temperature change during combustion (pre-period, firing, main-period, post-period).
3. Calculating the heat of combustion (H) using the relations:

$$
Q = C_{cal} \cdot \Delta T \qquad H = \frac{Q}{m}
$$

4. Performing uncertainty and error propagation:

$$
\frac{\Delta Q}{Q} = \sqrt{\left(\frac{\Delta C_{cal}}{C_{cal}}\right)^2 + \left(\frac{u_{\Delta T}}{\Delta T}\right)^2}
$$

$$
\frac{\Delta H}{H} = \sqrt{\left(\frac{\Delta Q}{Q}\right)^2 + \left(\frac{u_m}{m}\right)^2}
$$

where $u_{\Delta T}$ is the instrumental uncertainty of the temperature sensor and $u_m$ is the uncertainty of the analytical balance (not to be confused with $\Delta T$ itself).

## Standards Alignment

The calculation modules are structured to be extensible toward major international calorimetry standards, including:

* **Liquid Fuels:** ASTM D240, ASTM D4809
* **Solid Mineral Fuels & Coal:** ISO 1928, ASTM D5865
* **Biomass & Biofuels:** ISO 18125
* **Waste & Refuse-Derived Fuels:** ASTM D5468, ASTM E711

> **Note:** The current implementation calculates ΔT as a simple difference between the firing-point and maximum vessel temperatures. It does **not** yet apply the heat-leak / radiation correction (e.g. Regnault-Pfaundler method) required for full compliance with these standards. Pre-period and post-period drift data are parsed but not yet used in the correction. Full standards compliance is planned but not yet implemented.

## Tech Stack

* Python 3
* Pandas (data I/O and manipulation)
* Matplotlib (visualization)
* SciPy (optimization and curve fitting; used for future correction methods)

## Installation

```bash
git clone <repo-url>
cd bomb-calorimeter-analysis
pip install -r requirements.txt
```

## Usage

```bash
python processor.py
```

Place your instrument export file at `data/sample_data.csv`. The file should contain a metadata header (lines starting with `*`) followed by CSV sensor data:

```
* =========================================
* DEVICE / DATE / SAMPLE ID / SAMPLE MASS / TEST STANDARD
* =========================================
Time(s),Vessel_Temp(C),Jacket_Temp(C),Phase
0,25.1012,25.0001,Pre-Period
...
```

The `Phase` column must include a `Firing` row marking ignition, plus `Pre-Period`, `Main-Period`, and `Post-Period` phases.

## Results & Error Analysis

The thermodynamic calculations incorporate error propagation based on GUM (Guide to the Expression of Uncertainty in Measurement) principles. Accounting for the instrumental uncertainties of the temperature sensor (± 0.0005 °C), the analytical balance (± 0.0001 g), and the calibration constant, the specific heat of combustion for the sample dataset (benzoic acid standard) was determined as:

* **Final Result:** H = 26473.4 ± 37.2 J/g
* **Relative Error:** 0.141 %

A relative uncertainty below 0.2% reflects the precision of the simulated IC210 system data and the processing pipeline, given the current (uncorrected) ΔT calculation. Applying a heat-leak correction is expected to shift the result and is the next planned step.

## License

<!-- Add a license, e.g. MIT -->
