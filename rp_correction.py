"""
Regnault-Pfaundler Heat-Leak Correction
=========================================
Implements the classical Regnault-Pfaundler (R-P) correction for isoperibol
bomb calorimetry, per the method referenced in ASTM D240 / ISO 1928.

The vessel exchanges heat with its surroundings throughout the whole
experiment (not just during the main/reaction period). The Pre-Period and
Post-Period drift rates are used to estimate this continuous heat leak and
subtract it from the raw temperature rise, giving a corrected Delta T that
better represents the temperature rise attributable to combustion alone.

Formula
-------
    a  = firing (ignition) time
    T_a = vessel temperature at time a
    c  = time at which the Main-Period ends / Post-Period begins
    T_c = vessel temperature at time c
    b  = time at which the vessel reaches 60% of the total (T_c - T_a) rise
         (approximates the inflection point of the temperature-time curve)
    r1 = drift rate (deg/s) fitted over the Pre-Period (before firing)
    r2 = drift rate (deg/s) fitted over the Post-Period (after Main-Period)

    Delta T_corrected = (T_c - T_a) - r1*(b - a) - r2*(c - b)

This corrected Delta T should be used in Q = C_cal * Delta T instead of the
raw (T_final - T_initial) difference.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import linregress


def load_sensor_data(file_path):
    df = pd.read_csv(file_path, comment='*')
    df.columns = df.columns.str.strip()
    return df


def fit_drift_rate(df, phase_label, time_col='Time(s)', temp_col='Vessel_Temp(C)'):
    """Linear-fit the temperature drift rate (deg/s) over a given phase."""
    phase_df = df[df['Phase'] == phase_label]
    if len(phase_df) < 2:
        raise ValueError(f"Not enough points in phase '{phase_label}' to fit a rate.")
    result = linregress(phase_df[time_col], phase_df[temp_col])
    return result.slope, result.intercept, result.rvalue ** 2


def find_time_at_fraction(df, t_a, T_a, t_c, T_c, fraction=0.6,
                           time_col='Time(s)', temp_col='Vessel_Temp(C)'):
    """
    Find time 'b' at which the vessel temperature first reaches the given
    fraction of the total rise between firing (a) and end of main period (c).
    Uses linear interpolation between the bracketing samples.
    """
    target_T = T_a + fraction * (T_c - T_a)

    window = df[(df[time_col] >= t_a) & (df[time_col] <= t_c)].sort_values(time_col)
    times = window[time_col].to_numpy()
    temps = window[temp_col].to_numpy()

    if target_T <= temps[0]:
        return times[0]
    if target_T >= temps[-1]:
        return times[-1]

    idx = np.searchsorted(temps, target_T)
    t0, t1 = times[idx - 1], times[idx]
    T0, T1 = temps[idx - 1], temps[idx]
    # linear interpolation
    b = t0 + (target_T - T0) * (t1 - t0) / (T1 - T0)
    return b


def regnault_pfaundler_correction(df, time_col='Time(s)', temp_col='Vessel_Temp(C)'):
    """
    Computes the R-P corrected Delta T for a single calorimeter run.

    Returns a dict with the corrected Delta T, the raw (uncorrected) Delta T,
    fitted drift rates, and the key time/temperature markers, so the result
    is fully auditable.
    """
    # --- Key markers ---
    t_a = df[df['Phase'] == 'Firing'][time_col].iloc[0]
    T_a = df[df['Phase'] == 'Firing'][temp_col].iloc[0]

    post_period = df[df['Phase'] == 'Post-Period']
    if post_period.empty:
        raise ValueError("No Post-Period data found; cannot apply R-P correction.")
    t_c = post_period[time_col].iloc[0]
    T_c = df[df[time_col] <= t_c][temp_col].max()  # peak temp up to start of post-period

    # --- Drift rates from Pre- and Post-Period ---
    r1, _, r1_r2 = fit_drift_rate(df, 'Pre-Period', time_col, temp_col)
    r2, _, r2_r2 = fit_drift_rate(df, 'Post-Period', time_col, temp_col)

    # --- Inflection-point time (60% rise convention) ---
    t_b = find_time_at_fraction(df, t_a, T_a, t_c, T_c, fraction=0.6,
                                 time_col=time_col, temp_col=temp_col)

    # --- Raw and corrected Delta T ---
    delta_T_raw = T_c - T_a
    leak_before = r1 * (t_b - t_a)
    leak_after = r2 * (t_c - t_b)
    delta_T_corrected = delta_T_raw - leak_before - leak_after

    return {
        "t_a": t_a, "T_a": T_a,
        "t_b": t_b,
        "t_c": t_c, "T_c": T_c,
        "r1_per_s": r1, "r1_fit_r2": r1_r2,
        "r2_per_s": r2, "r2_fit_r2": r2_r2,
        "delta_T_raw": delta_T_raw,
        "leak_correction_before_b": leak_before,
        "leak_correction_after_b": leak_after,
        "delta_T_corrected": delta_T_corrected,
    }


def main():
    file_path = os.path.join("data", "sample_data.csv")
    C_cal = 10815.0  # J/deg
    sample_mass = 0.9845  # g

    df = load_sensor_data(file_path)
    result = regnault_pfaundler_correction(df)

    Q_raw = C_cal * result["delta_T_raw"]
    Q_corr = C_cal * result["delta_T_corrected"]
    H_raw = Q_raw / sample_mass
    H_corr = Q_corr / sample_mass

    print("=" * 55)
    print("   REGNAULT-PFAUNDLER HEAT-LEAK CORRECTION REPORT")
    print("=" * 55)
    print(f"Firing time (a)        : t = {result['t_a']:.1f} s, T = {result['T_a']:.4f} C")
    print(f"Inflection time (b)    : t = {result['t_b']:.2f} s  (60% rise point)")
    print(f"Main-Period end (c)    : t = {result['t_c']:.1f} s, T = {result['T_c']:.4f} C")
    print("-" * 55)
    print(f"Pre-Period drift  (r1) : {result['r1_per_s']:.6f} C/s  (R^2={result['r1_fit_r2']:.4f})")
    print(f"Post-Period drift (r2) : {result['r2_per_s']:.6f} C/s  (R^2={result['r2_fit_r2']:.4f})")
    print("-" * 55)
    print(f"Raw Delta T            : {result['delta_T_raw']:.4f} C")
    print(f"  - leak correction (a-b): {result['leak_correction_before_b']:+.4f} C")
    print(f"  - leak correction (b-c): {result['leak_correction_after_b']:+.4f} C")
    print(f"Corrected Delta T       : {result['delta_T_corrected']:.4f} C")
    print("=" * 55)
    print(f"H (raw, uncorrected)    : {H_raw:.1f} J/g")
    print(f"H (R-P corrected)       : {H_corr:.1f} J/g")
    print(f"Difference              : {H_corr - H_raw:+.1f} J/g "
          f"({100 * (H_corr - H_raw) / H_raw:+.3f} %)")
    print("=" * 55)


if __name__ == "__main__":
    main()
