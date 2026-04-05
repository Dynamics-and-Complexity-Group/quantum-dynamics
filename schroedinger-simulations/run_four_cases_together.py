# Last updated on March 8, 2026

import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags
from scipy.integrate import simps
import time
from multiprocessing import Process
import matplotlib.pyplot as plt

# ============================================================
# -------------------- Physical Quantities -------------------
# ============================================================

def energy(psi, dx, V):
    d2psi = np.zeros_like(psi)
    d2psi[1:-1] = (psi[2:] - 2*psi[1:-1] + psi[:-2]) / dx**2
    kinetic = -0.5 * np.sum(np.conj(psi) * d2psi) * dx
    potential = np.sum(V * np.abs(psi)**2) * dx
    return kinetic + potential


def compute_moments(psi, x, dx):
    prob = np.abs(psi)**2
    mean = np.sum(x * prob) * dx
    var = np.sum((x - mean)**2 * prob) * dx
    std = np.sqrt(var)

    if std == 0:
        return mean, var, 0.0, 0.0

    skew = np.sum((x - mean)**3 * prob) * dx / std**3
    kurt = np.sum((x - mean)**4 * prob) * dx / std**4 - 3  # excess kurtosis

    return mean, var, skew, kurt

def plot_moments(times, mean, variance, skewness, kurtosis,
                 filename="moments_plot.png", dpi=150):
    """
    Plot mean position, variance, skewness, and excess kurtosis vs time.
    """

    fig, axes = plt.subplots(4, 1, figsize=(8, 12), sharex=True)

    # --- Mean ---
    axes[0].plot(times, mean, 'b-', linewidth=2)
    axes[0].set_ylabel("Mean ⟨x⟩")
    axes[0].axhline(0, color="gray", linestyle="--", alpha=0.5)

    # --- Variance ---
    axes[1].plot(times, variance, 'r-', linewidth=2)
    axes[1].set_ylabel("Variance σ²")

    # --- Skewness ---
    axes[2].plot(times, skewness, 'g-', linewidth=2)
    axes[2].axhline(0, color="gray", linestyle="--", alpha=0.5)
    axes[2].set_ylabel("Skewness")

    # --- Excess Kurtosis ---
    axes[3].plot(times, kurtosis, 'm-', linewidth=2)
    axes[3].axhline(0, color="gray", linestyle="--", alpha=0.5)
    axes[3].set_ylabel("Excess Kurtosis")
    axes[3].set_xlabel("Time")

    plt.tight_layout()
    plt.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    print(f"Plot saved to '{filename}'")

# ============================================================
# -------------------- Main Solver ----------------------------
# ============================================================

def main(well="left", x_in=1.0, var_val=0.4, E_target=12.0):

    print(f"\nRunning {well.upper()} well simulation")
    
    # --- Units ---
    hbar, m = 1.0, 1.0

    # --- Grid ---
    x_min, x_max, N = -100, 100, 10**5
    dx = (x_max - x_min) / (N - 1)
    x = np.linspace(x_min, x_max, N)

    dt = 0.01
    num_steps = 15 * 1024
    times = np.arange(num_steps + 1) * dt

    # --- Potential ---
    a, b, c = 10.0, 4.0, 0.35
    V = a/2*x**2 - b/3*x**3 + c/4*x**4

    # --- Well-dependent settings ---
    if well == "left":
        momentum_sign = +1
        prefix = "left_well"
    elif well == "right":
        momentum_sign = -1
        prefix = "right_well"
    else:
        raise ValueError("well must be 'left' or 'right'")

    sigma = np.sqrt(var_val)
    x0 = x_in

    V0 = V[np.argmin(np.abs(x - x0))]
    kinetic_needed = E_target - V0
    k = np.sqrt(2 * kinetic_needed)

    psi = (1 / (2*np.pi*sigma**2))**0.25 \
        * np.exp(-(x - x0)**2 / (4*sigma**2)) \
        * np.exp(1j * momentum_sign * k * (x - x0))

    psi /= np.sqrt(simps(np.abs(psi)**2, x))

    E_initial = energy(psi, dx, V)
    print(f"Computed initial energy: {E_initial.real:.6f}")


    # --- Allocate moment arrays ---
    means = np.zeros(num_steps + 1)
    variances = np.zeros(num_steps + 1)
    skewnesses = np.zeros(num_steps + 1)
    kurtoses = np.zeros(num_steps + 1)

    means[0], variances[0], skewnesses[0], kurtoses[0] = compute_moments(psi, x, dx)

    # --- Crank–Nicolson matrices ---
    alpha = 1j * dt * hbar / (4 * m * dx**2)

    A_diag = 1 + 2*alpha + (1j * dt * V)/(2*hbar)
    B_diag = 1 - 2*alpha - (1j * dt * V)/(2*hbar)
    A_off = -alpha * np.ones(N - 1, dtype=complex)
    B_off = alpha * np.ones(N - 1, dtype=complex)

    A = diags([A_off, A_diag, A_off], [-1, 0, 1], format='csr')
    B = diags([B_off, B_diag, B_off], [-1, 0, 1], format='csr')

    print("Starting time evolution...")

    for step in range(num_steps):

        psi = spsolve(A, B.dot(psi))
        psi[0] = psi[-1] = 0

        means[step+1], variances[step+1], \
        skewnesses[step+1], kurtoses[step+1] = compute_moments(psi, x, dx)

        if (step + 1) % 100 == 0:
            current_norm = np.sum(np.abs(psi)**2) * dx
            current_energy = energy(psi, dx, V)
            print(f"{well.upper()} | Step {step+1}: "
                  f"Norm Dev = {current_norm - 1:.2e}, "
                  f"Energy Dev = {(current_energy - E_initial).real:.2e}")

    final_energy = energy(psi, dx, V)
    energy_str = f"{final_energy.real:.4f}".replace('.', '_')

    # --- Save data (moments only) ---
    np.savez(f"asymetric_double_well_{prefix}_{energy_str}.npz",
             times=times,
             mean=means,
             variance=variances,
             skewness=skewnesses,
             kurtosis=kurtoses,
             E_initial=E_initial,
             E_final=final_energy,
             x_in=x_in,
             var_val=var_val,
             E_target=E_target)

    print(f"{well.upper()} well complete.\n")

    # --- Plot moments ---
    plot_moments(times,means,variances,skewnesses, kurtoses,
        filename=f"moments_{prefix}_{energy_str}.png"
    )

# ============================================================
# -------------------- Parallel Execution ---------------------
# ============================================================

if __name__ == "__main__":

    start_time = time.time()

    processes = []

    # LEFT well: below & above threshold
    processes.append(Process(target=main,
        kwargs=dict(well="left", x_in=1.0, var_val=0.45, E_target=12.66)))
    processes.append(Process(target=main,
        kwargs=dict(well="left", x_in=0.5, var_val=0.45, E_target=17.89)))

    # RIGHT well: below & above threshold
    processes.append(Process(target=main,
        kwargs=dict(well="right", x_in=5.5, var_val=0.1, E_target=12.55)))
    processes.append(Process(target=main,
        kwargs=dict(well="right", x_in=5.5, var_val=0.1, E_target=18.5)))

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    elapsed = (time.time() - start_time) / 60
    print(f"\nTotal execution time: {elapsed:.4f} minutes")