"""
This script numerically solves the time-dependent Schrödinger equation
for a particle in an asymmetric double-well potential using the
Crank–Nicolson method. The simulation visualizes the quantum tunneling
process from the left to the right well.

Output: Animation saved as MP4 file and plots of mean position and variance over time.
(saving data fileS disabled for now)

Last updated: October 24, 2025.
"""

import numpy as np
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import time
from scipy.sparse import diags
from scipy.integrate import simps
import pandas as pd
from matplotlib.animation import FuncAnimation  # Import FuncAnimation for animation

def energy(psi, dx, V, hbar=1.0, m=1.0):
    """
    More consistent energy calculation using finite differences
    """
    # Kinetic energy (using second-order finite differences)
    d2psi = np.zeros_like(psi)
    d2psi[1:-1] = (psi[2:] - 2*psi[1:-1] + psi[:-2])/dx**2
    kinetic = -0.5 * np.sum(np.conj(psi) * d2psi) * dx
    potential = np.sum(V * np.abs(psi)**2) * dx
    return kinetic + potential

def mean_position(psi, dx, x):
    """
    Calculates the expectation value of position (mean position) for a wave function.
    
    Parameters:
    psi : numpy array
        The wave function values on a discrete grid.
    dx : float
        The grid spacing.
    x : numpy array
        The position grid points corresponding to psi.
    
    Returns:
    float
        The expectation value of position ⟨x⟩.
    """
    probability_density = np.abs(psi)**2
    mean_pos = np.sum(x * probability_density) * dx
    return mean_pos

def variance_position(psi, dx, x):
    """
    Calculates the variance of position for a wave function.
    
    Parameters:
    psi : numpy array
        The wave function values on a discrete grid.
    dx : float
        The grid spacing.
    x : numpy array
        The position grid points corresponding to psi.
    
    Returns:
    float
        The variance of position ⟨(x - ⟨x⟩)²⟩.
    """
    # First calculate the mean position
    probability_density = np.abs(psi)**2
    mean_pos = np.sum(x * probability_density) * dx
    
    # Calculate the variance: ⟨(x - ⟨x⟩)²⟩
    variance = np.sum((x - mean_pos)**2 * probability_density) * dx
    return variance

def skewness_position(psi, dx, x):
    """
    Calculates the skewness of the position distribution for a wave function.
    
    Parameters:
    psi : numpy array
        The wave function values on a discrete grid.
    dx : float
        The grid spacing.
    x : numpy array
        The position grid points corresponding to psi.
    
    Returns:
    float
        The skewness of the position distribution:
        ⟨[(x - ⟨x⟩)³]⟩ / (⟨[(x - ⟨x⟩)²]⟩)^(3/2)
    """
    probability_density = np.abs(psi)**2
    mean_pos = np.sum(x * probability_density) * dx
    variance = np.sum((x - mean_pos)**2 * probability_density) * dx
    std_dev = np.sqrt(variance)
    
    # Calculate the third central moment
    third_moment = np.sum((x - mean_pos)**3 * probability_density) * dx
    
    # Avoid division by zero if variance is zero
    if std_dev == 0:
        return 0.0
    
    skewness = third_moment / (std_dev**3)
    return skewness


def create_animation(psi_data, x, V, dx, dt, filename='double_well_animation.mp4', fps=50, dpi=150):
    """
    Creates and saves an animation of the probability density, potential, and mean position.

    Args:
        psi_data (np.ndarray): Array containing psi at each time step (num_steps + 1, N).
        x (np.ndarray): Spatial grid.
        V (np.ndarray): Potential array.
        dx (float): Spatial grid spacing.
        dt (float): Time step.
        filename (str): Name of the output video file.
        fps (int): Frames per second for the animation.
        dpi (int): Dots per inch for video resolution.
    """
    print(f"\nCreating animation with mean position. This may take some time depending on num_steps and fps...")
    fig, ax = plt.subplots(figsize=(12, 7))

    # Initialize plot lines with the first frame's data
    line_psi_sq, = ax.plot(x, np.abs(psi_data[0,:])**2, label='$|\psi|^2$')
    line_potential, = ax.plot(x, V*0.05, 'g:', label='Potential $V(x)$')
    mean_pos = mean_position(psi_data[0, :], dx, x)
    point_mean, = ax.plot(mean_pos, 0, 'ro', label='Mean Position ⟨x⟩')
    ax.axhline(0, color='k', linestyle='--', alpha=0.5)  # Zero line for reference
    ax.axvline(x=0.0, color="gray", linestyle="--")
    ax.axvline(x=3.693980625181293, color="gray", alpha=0.5, linestyle="--")
    ax.axvline(x=7.7345908033901365, color="gray", linestyle="--")
    ax.set_xlabel("Position (x)")
    ax.set_ylabel("Probability Density")
    ax.legend(loc='upper center')
    ax.set_xlim(-4.0, 10.0)
    ax.set_ylim(-1.0, 4.0)

    def update(frame):
        """Update function for each frame of the animation."""
        current_psi_sq = np.abs(psi_data[frame,:])**2
        line_psi_sq.set_ydata(current_psi_sq)
        mean_pos = mean_position(psi_data[frame, :], dx, x)
        point_mean.set_xdata([mean_pos])
        ax.set_title(f"(t={frame * dt:.3f})")
        return line_psi_sq, point_mean  # Return all artists that were modified

    anim = FuncAnimation(fig, update, frames=range(0, len(psi_data)), interval=1000/fps, blit=False)

    try:
        # Requires ffmpeg installed and in your system's PATH
        anim.save(filename, writer='ffmpeg', dpi=dpi)
        print(f"Animation saved successfully as '{filename}'.")
    except Exception as e:
        print(f"Error saving animation: {e}")
        print("Please ensure 'ffmpeg' is installed on your system and accessible in your PATH.")
        print("You can download ffmpeg from: https://ffmpeg.org/download.html")
        print("For Windows, download the zip, extract, and add the /bin folder to your system's PATH environmental variable.")

    plt.close(fig)  # Close the animation figure to prevent it from lingering after saving

def plot_mean_position_and_variance_and_skewness(psi_data, x, dx, dt, filename='mean_position_variance_skewness_plot.png', dpi=150):
    """
    Plots and saves the mean position and variance as functions of time.

    Args:
        psi_data (np.ndarray): Array containing psi at each time step (num_steps + 1, N).
        x (np.ndarray): Spatial grid.
        dx (float): Spatial grid spacing.
        dt (float): Time step.
        filename (str): Name of the output image file.
        dpi (int): Dots per inch for image resolution.
    """
    print(f"\nCreating mean position and variance plot...")
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Compute times, mean positions, and variances
    times = np.arange(len(psi_data)) * dt
    mean_positions = np.array([mean_position(psi_data[i, :], dx, x) for i in range(len(psi_data))])
    variances = np.array([variance_position(psi_data[i, :], dx, x) for i in range(len(psi_data))])
    skewnesses = np.array([skewness_position(psi_data[i, :], dx, x) for i in range(len(psi_data))])

    
    # Create DataFrame and save as CSV
    df = pd.DataFrame({
        'time': times, 
        'mean_position': mean_positions,
        'variance': variances,
        'skewness':skewnesses
    })

    # save data as csv
    #df.to_csv('time_position_variance_data_right_well_Energy_23_442.csv', index=False)
    # df.to_csv('time_position_variance_data_right_well_Energy_21_9421.csv', index=False)
    
   # Create 3x1 subplot layout
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    # --- Plot mean position vs time (top subplot) ---
    ax1.plot(times, mean_positions, 'b-', label='Mean Position ⟨x⟩', linewidth=2)
    ax1.set_ylabel("Mean Position ⟨x⟩")
    ax1.set_title("Mean Position, Variance, and Skewness over Time")
    ax1.legend(loc='upper right')
    ax1.axhline(y=0.0, color="gray", linestyle="--")
    ax1.axhline(y=3.6939806251812937, color="gray", alpha=0.5, linestyle="--")
    ax1.axhline(y=7.7345908033901365, color="gray", linestyle="--")
    ax1.set_ylim([-0.30, 10.2])

    # --- Plot variance vs time (middle subplot) ---
    ax2.plot(times, variances, 'r-', label='Variance σ²', linewidth=2)
    ax2.set_ylabel("Variance σ²")
    ax2.legend(loc='upper right')
    # ax2.grid(True, alpha=0.3)

    # --- Plot skewness vs time (bottom subplot) ---
    ax3.plot(times, skewnesses, 'g-', label='Skewness', linewidth=2)
    ax3.set_xlabel(r"Time $t$")
    ax3.set_ylabel("Skewness")
    ax3.legend(loc='upper right')
    ax3.axhline(y=0.0, color="gray", linestyle="--", alpha=0.6)
    # ax3.grid(True, alpha=0.3)

    # --- Adjust layout and save ---
    plt.tight_layout()
    plt.savefig("mean_variance_skewness_vs_time.png", dpi=300)


        
    try:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"Mean position and variance plot saved successfully as '{filename}'.")
    except Exception as e:
        print(f"Error saving plot: {e}")
        
    plt.close(fig)  # Close the figure to prevent it from lingering
        
    # Also save the data as NPZ file for numerical analysis
    # --- Save Data ---
    #data_filename = filename.replace(".png", ".npz")

    #np.savez(data_filename, 
            #times=times, 
            #mean_positions=mean_positions, 
            #variances=variances)
    #print("Data also saved.")


def main(x_in = 0.5, var_val = 0.5, E_target=12):

    print("computing for E_target = ", E_target)
    # --- Parameters for a Symmetric Double Well ---
    hbar, m = 1.0, 1.0  # Natural units
    #omega, lambda_ = np.sqrt(20), 20.0

    # Spatial domain, N grid points
    x_min, x_max, N = -100, 100, 10**5  # Larger domain, high N for resolution
    dx = (x_max - x_min) / (N - 1)  # Spatial step

    # Time step
    dt = 0.01
    num_steps = 15*1024

    alpha = 1j * dt * hbar / (4 * m * dx**2)  # Crank-Nicolson coefficient

    # Spatial grid
    x = np.linspace(x_min, x_max, N)

    # Define asymmetric double well potential
    a = 10.0
    b = 4.0
    c = 0.35 
    V = a/2*x**2 - b/3*x**3 + c/4*x**4

    # --- Initial wave function: Gaussian wave packet ---
    # Parameters for quantum stability study

    x0 = x_in

    # Left well minima at x = 0

    # Barrier at x = 3.693980625181293
    # Barrier height = 17.3117

    # Righ well minima at 7.7346
    # Right well depth 4.6780

    # Total height of the barrier 17.3117 + 4.6780 = 21.9897

    # threshold from fixed point analysis 
    # E_min = 8.531237506247502
    # For instability, E > 10.60500640208603
    # E_threshold = 10.60500640208603 + 4.6780 = 15.28300640208603

    # For tunneling: traget 17.678 (above 15.28)
    # No  tunneling: traget 13.678 (below 15.28)

    E_target =  E_target

    # Choose an appropriate sigma (width of Gaussian)
    sigma = np.sqrt(var_val)  

    # Calculate the potential energy at the center of the wave packet
    V0 = V[np.argmin(np.abs(x - x0))]  # Potential at x0

    #print(V0)
    # Calculate the kinetic energy needed to reach E_target
    # Total energy = kinetic + potential
    # So kinetic = E_target - potential
    kinetic_needed = E_target - V0

    # Calculate k from kinetic energy: KE = (ħ²k²)/(2m)
    # With ħ=1, m=1: KE = k²/2 → k = sqrt(2*KE)
    k = np.sqrt(2 * kinetic_needed)


    print("wave vector k =", k)

    psi = (1 / (2.0*np.pi * sigma**2))**0.25 * np.exp(-(x - x0)**2 / (4 * sigma**2)) * np.exp(1j * (k) * (x - x0))
    psi /= np.sqrt(simps(np.abs(psi)**2, x))  # Normalize

    # Calculate initial energy
    E_initial = energy(psi, dx, V, hbar, m)
    print(f"Calculated initial energy (numerical): {E_initial.real:.4f} (imag: {E_initial.imag:.2e})")

    # Additional diagnostics for your stability study
    mean_x = mean_position(psi, dx, x)
    print(f"Initial mean position: {mean_x:.4f}")
    print(f"Distance to barrier: {abs(mean_x - 3.60):.4f}")

    # Wave packet width analysis
    psi_sq = np.abs(psi)**2
    var_x = np.sum((x - mean_x)**2 * psi_sq) * dx
    print(f"Initial wave packet width (σ): {np.sqrt(var_x):.4f}")

    # Store wave function at each time step
    psi_data = np.zeros((num_steps + 1, N), dtype=complex)
    psi_data[0, :] = psi.copy()

    # Diagonal and off-diagonal elements for sparse matrices
    A_diag = 1 + 2 * alpha + (1j * dt * V) / (2 * hbar)
    A_off = -alpha * np.ones(N - 1, dtype=complex)

    B_diag = 1 - 2 * alpha - (1j * dt * V) / (2 * hbar)
    B_off = alpha * np.ones(N - 1, dtype=complex)

    # Construct tridiagonal sparse matrices
    A = diags([A_off, A_diag, A_off], [-1, 0, 1], shape=(N, N), dtype=complex, format='csr')
    B = diags([B_off, B_diag, B_off], [-1, 0, 1], shape=(N, N), dtype=complex, format='csr')

    # --- Main Simulation Loop ---
    print("\nStarting time evolution...")
    for step in range(num_steps):
        psi = spsolve(A, B.dot(psi))
        psi[0] = psi[-1] = 0  # Enforce boundary conditions (zero at boundaries)

        psi_data[step+1, :] = psi.copy()  # Store psi for animation/analysis

        # --- Monitoring Output ---
        current_norm = np.sum(np.abs(psi)**2) * dx
        if abs(current_norm - 1) > 1e-10:
            print(f"Step {step+1}: Norm deviation: {current_norm - 1:.2e}")

        if (step + 1) % 100 == 0:  # Check energy less frequently for long simulations
            current_energy = energy(psi, dx, V, hbar, m)
            energy_deviation = current_energy - E_initial
            print(f"Step {step+1}: Norm = {current_norm:.6f} (Dev: {current_norm - 1:.2e}), Energy = {current_energy.real:.4f} (Dev: {energy_deviation.real:.2e})")

    print("\nTime evolution complete.")

    # --- Final Checks ---
    final_norm = np.sum(np.abs(psi)**2) * dx
    final_energy = energy(psi, dx, V, hbar, m)
    final_energy_deviation = final_energy - E_initial

    print(f"Final Norm: {final_norm:.6f} (Deviation: {final_norm - 1:.2e})")
    print(f"Final Energy: {final_energy.real:.4f} (Deviation: {final_energy_deviation.real:.2e})")

    # --- Save Data ---
    # Convert the energy to string and replace '.' with '_'
    energy_str = f"{final_energy.real:.4f}".replace('.', '_')

    #data_filename = f"asymetric_double_well_simulation_data_right_well_{energy_str}.npz"
    #np.savez(data_filename, psi_data=psi_data, x=x, V=V, dt=dt, E_initial=E_initial)
    #print(f"\nSimulation data saved to '{data_filename}'.")

    # --- Create and Save Animation ---

    animation_filename =  f"asymetric_double_well_simulation_data_right_well_{energy_str}.mp4"
    create_animation(psi_data, x, V, dx, dt, filename=animation_filename, fps=30, dpi=150)

    # --- Plot and Save Mean Position, Variance, and Skewness ---
    plot_mean_position_and_variance_and_skewness(psi_data, x, dx, dt, filename=f'asymetric_double_well_simulation_data_left_well_{energy_str}.png', dpi=150)

#################################################################################################################################################
#################################################################################################################################################


if __name__ == "__main__":
    start_time = time.time()  # record start time
    print("Expecte execution time: ~30 minutes for 15k steps with 100k grid points.\n")
    main(x_in = 0.5, var_val = 0.45, E_target = 11.943) # Gives E = 13.6869 =4.68 + 9.0069 # No tunneling
    main(x_in = 0.5, var_val = 0.45, E_target = 17.89) # Gives E = 19.6338 = 4.68 + 14.9538 # Tunneling 
    end_time = time.time()    # record end time
    elapsed_minutes = (end_time - start_time) / 60
    print(f"Execution time: {elapsed_minutes:.4f} minutes")

