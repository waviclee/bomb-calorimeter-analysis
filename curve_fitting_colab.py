import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress

plt.style.use("ggplot")

# =====================================================================
# 1. Ensure sample data exists (creates it if running in a fresh Colab)
# =====================================================================
os.makedirs("data", exist_ok=True)
DATA_PATH = "data/sample_data.csv"

if not os.path.exists(DATA_PATH):
    csv_content = """* =========================================
* DEBYE TECHNIC THERMAL ANALYSIS SYSTEMS
* Device: IC210 Isoperibol Calorimeter
* Date: 2026-07-16
* Sample ID: Benzoic_Acid_Std
* Sample Mass (g): 0.9845
* Test Standard: ASTM D240
* =========================================
Time(s),Vessel_Temp(C),Jacket_Temp(C),Phase
0,25.1012,25.0001,Pre-Period
10,25.1021,25.0001,Pre-Period
20,25.1038,25.0002,Pre-Period
30,25.1051,25.0001,Pre-Period
40,25.1065,25.0001,Pre-Period
50,25.1082,25.0002,Pre-Period
60,25.1105,25.0001,Firing
70,25.4503,25.0001,Main-Period
80,26.1208,25.0002,Main-Period
90,26.8504,25.0001,Main-Period
100,27.2109,25.0001,Main-Period
110,27.4502,25.0002,Main-Period
120,27.5204,25.0001,Main-Period
130,27.5158,25.0002,Main-Period
140,27.5101,25.0001,Main-Period
150,27.5024,25.0001,Main-Period
160,27.4957,25.0001,Main-Period
170,27.4882,25.0002,Main-Period
180,27.4815,25.0001,Main-Period
190,27.4741,25.0001,Main-Period
200,27.4658,25.0002,Main-Period
210,27.4582,25.0001,Post-Period
220,27.4501,25.0001,Post-Period
230,27.4435,25.0002,Post-Period
240,27.4362,25.0001,Post-Period
250,27.4298,25.0001,Post-Period
260,27.4221,25.0002,Post-Period
270,27.4154,25.0001,Post-Period
280,27.4082,25.0001,Post-Period
290,27.4015,25.0002,Post-Period
300,27.3951,25.0001,Post-Period"""
    with open(DATA_PATH, "w", encoding="utf-8") as f:
        f.write(csv_content)
    print(f"'{DATA_PATH}' not found, created a sample file.")


def load_sensor_data(file_path):
    df = pd.read_csv(file_path, comment='*')
    df.columns = df.columns.str.strip()
    return df


# =====================================================================
# 2. Exponential cooling model: T(t) = T_inf + A * exp(-k * (t - t0))
# =====================================================================
def exponential_model(t, T_inf, A, k, t0):
    return T_inf + A * np.exp(-k * (t - t0))


def fit_exponential(df, phase_label, time_col='Time(s)', temp_col='Vessel_Temp(C)'):
    """
    Fits the exponential cooling model to one phase's data.
    Returns fitted parameters, their 1-sigma uncertainties, and fit R^2
    (also returns the linear-fit R^2 for comparison).
    """
    phase_df = df[df['Phase'] == phase_label].sort_values(time_col)
    t = phase_df[time_col].to_numpy(dtype=float)
    T = phase_df[temp_col].to_numpy(dtype=float)

    if len(t) < 4:
        raise ValueError(f"Need at least 4 points to fit phase '{phase_label}'.")

    t0 = t[0]
    T_inf0 = T[-1]
    A0 = T[0] - T[-1]
    k0 = 1e-3

    popt, pcov = curve_fit(
        lambda tt, T_inf, A, k: exponential_model(tt, T_inf, A, k, t0),
        t, T, p0=[T_inf0, A0, k0], maxfev=10000
    )
    T_inf, A, k = popt
    perr = np.sqrt(np.diag(pcov))

    T_pred = exponential_model(t, T_inf, A, k, t0)
    ss_res = np.sum((T - T_pred) ** 2)
    ss_tot = np.sum((T - np.mean(T)) ** 2)
    r2_exp = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    lin = linregress(t, T)
    r2_lin = lin.rvalue ** 2

    return {
        "t": t, "T": T, "t0": t0,
        "T_inf": T_inf, "A": A, "k": k,
        "T_inf_err": perr[0], "A_err": perr[1], "k_err": perr[2],
        "r2_exponential": r2_exp, "r2_linear": r2_lin,
        "linear_slope": lin.slope, "linear_intercept": lin.intercept,
    }


def instantaneous_rate(fit, t):
    """dT/dt of the fitted exponential curve at time t (deg C / s)."""
    return -fit["k"] * fit["A"] * np.exp(-fit["k"] * (t - fit["t0"]))


# =====================================================================
# 3. Run the fit for Pre-Period and Post-Period, print a report
# =====================================================================
df = load_sensor_data(DATA_PATH)
fits = {}

print("=" * 60)
print("   EXPONENTIAL COOLING FIT vs LINEAR FIT COMPARISON")
print("=" * 60)

for phase in ["Pre-Period", "Post-Period"]:
    fit = fit_exponential(df, phase)
    fits[phase] = fit

    print(f"\n[{phase}]")
    print(f"  T_inf = {fit['T_inf']:.4f} +/- {fit['T_inf_err']:.4f} C")
    print(f"  A     = {fit['A']:.4f} +/- {fit['A_err']:.4f} C")
    print(f"  k     = {fit['k']:.6f} +/- {fit['k_err']:.6f} 1/s")
    print(f"  R^2 (exponential) = {fit['r2_exponential']:.6f}")
    print(f"  R^2 (linear)      = {fit['r2_linear']:.6f}")

    if fit["r2_exponential"] - fit["r2_linear"] > 1e-4:
        print("  -> Exponential model fits meaningfully better than linear.")
    else:
        print("  -> Linear approximation is adequate for this phase.")

print("\n" + "=" * 60)

# =====================================================================
# 4. Plot: raw data + fitted curves (exponential vs linear) side by side
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

for ax, phase in zip(axes, ["Pre-Period", "Post-Period"]):
    fit = fits[phase]
    t, T = fit["t"], fit["T"]
    t_dense = np.linspace(t.min(), t.max(), 200)

    T_exp = exponential_model(t_dense, fit["T_inf"], fit["A"], fit["k"], fit["t0"])
    T_lin = fit["linear_slope"] * t_dense + fit["linear_intercept"]

    ax.scatter(t, T, color='black', s=30, zorder=3, label='Raw data')
    ax.plot(t_dense, T_exp, color='firebrick', linewidth=2.2,
            label=f'Exponential fit (R²={fit["r2_exponential"]:.5f})')
    ax.plot(t_dense, T_lin, color='steelblue', linestyle='--', linewidth=2,
            label=f'Linear fit (R²={fit["r2_linear"]:.5f})')

    ax.set_title(phase, fontsize=12, fontweight='bold')
    ax.set_xlabel('Time (s)', fontweight='bold')
    ax.set_ylabel('Vessel Temperature (°C)', fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.6)

fig.suptitle('Exponential vs Linear Fit: Pre-Period and Post-Period Drift',
             fontsize=14, fontweight='bold')
fig.tight_layout()

os.makedirs("results/figures", exist_ok=True)
fig.savefig("results/figures/curve_fitting_comparison.png", dpi=300)
plt.show()

# =====================================================================
# 5. Example: instantaneous rate at the start/end of each phase
# =====================================================================
print("\nInstantaneous drift rate examples:")
for phase in ["Pre-Period", "Post-Period"]:
    fit = fits[phase]
    t_first, t_last = fit["t"][0], fit["t"][-1]
    r_first = instantaneous_rate(fit, t_first)
    r_last = instantaneous_rate(fit, t_last)
    print(f"  [{phase}] dT/dt at t={t_first:.0f}s: {r_first:.6f} C/s | "
          f"at t={t_last:.0f}s: {r_last:.6f} C/s | "
          f"(constant linear rate: {fit['linear_slope']:.6f} C/s)")
