import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PARAMETERS & REACTION KINETICS
# ==========================================
Lambda = 1.0          # population influx rate
alpha = 0.0002        # base environment-to-human transmission rate
beta = 0.0001         # base human to human transmission
mu = 0.0001           # natural death rate
omega = 0.01          # disease induced death rate
gamma = 0.2           # disease recovery rate (1/5 days from paper)
delta = 0.033         # removal rate of pathogen (1/30 days from paper)
zi = 0.3              # base shedding rate of infected hosts
g = 0.05              # intrinsic growth rate of the waterborne pathogen
K = 100000.0          # carrying capacity
d = 0.1               # pathogen diffusion rate (D4 in paper)
v = 0.5               # pathogen advection rate (river flow speed)
c_bound = 2.0         # relative bacterial loss rate at downstream (from paper pg 16)

# Diffusion rates (D1, D2, D3 from paper Page 15)
c1 = 1.0
c2 = 0.5
c3 = 0.8

def Reaction1(S, I, B, alpha_eff, beta_eff):
    return Lambda - alpha_eff*S*B - beta_eff*S*I - mu*S

def Reaction2(S, I, B, alpha_eff, beta_eff):
    return alpha_eff*S*B + beta_eff*S*I - (mu + omega + gamma)*I

def Reaction3(I, R):
    return gamma*I - mu*R

def Reaction4(B, I, zi_eff):
    return g*B * (1 - B/K) + zi_eff*I - delta*B


# ==========================================
# 2. THE NUMERICAL METHOD (SOLVER)
# ==========================================
def FiniteDiffSolver(L, T, nx, save_history=True, num_snapshots=100, fixed_dt=None):
    """
    Solves the 4-equation SIR-B model using FTCS finite difference.
    If fixed_dt is provided, it uses that exact time step (useful for accuracy checks).
    """
    dx = L / (nx - 1)
    
    # Time step calculation
    if fixed_dt is not None:
        nt = int(T / fixed_dt) + 1
        dt = T / nt # Recalculate slightly to hit T exactly
    else:
        # Calculate the maximum allowed dt based on the CFL condition for diffusion
        max_diffusion = max(c1, c2, c3, d)
        dt_max = 0.4 * (dx**2) / max_diffusion 
        nt = int(T / dt_max) + 1
        dt = T / nt

    # Spatial grid for initial conditions and spatial heterogeneity
    x = np.linspace(0, L, nx)
    H_x = 0.5 + 0.25 * np.cos(2 * np.pi * x) # Spatial Heterogeneity

    # Initial Conditions (Periodic)
    S = 10000 - 500 * np.cos(2 * x)
    I = 1 - np.cos(2 * x)
    R = np.zeros(nx)
    B = 0.5 - 0.3 * np.cos(2 * x)
    
    dS, dI, dR, dB = np.zeros(nx), np.zeros(nx), np.zeros(nx), np.zeros(nx)

    # History Tracking (Optional)
    if save_history:
        save_interval = max(1, nt // num_snapshots)
        S_hist, I_hist, R_hist, B_hist, t_hist = [], [], [], [], []

    # Time Loop
    for n in range(nt):
        current_time = n * dt
        
        if save_history and n % save_interval == 0:
            S_hist.append(S.copy()); I_hist.append(I.copy())
            R_hist.append(R.copy()); B_hist.append(B.copy())
            t_hist.append(current_time)

        # Seasonality Function
        T_t = 0.5 + 0.25 * np.sin(2 * np.pi * current_time / 12.0)
        
        # Effective rates
        alpha_eff = alpha * H_x * T_t
        beta_eff = beta * H_x * T_t
        zi_eff = zi * H_x * T_t

        # 2nd-Order Boundary Conditions (3-point one-sided differences)
        # Neumann BCs for human populations (dU/dx = 0)
        S[0] = (4*S[1] - S[2]) / 3.0
        S[-1] = (4*S[-2] - S[-3]) / 3.0
        
        I[0] = (4*I[1] - I[2]) / 3.0
        I[-1] = (4*I[-2] - I[-3]) / 3.0
        
        R[0] = (4*R[1] - R[2]) / 3.0
        R[-1] = (4*R[-2] - R[-3]) / 3.0
        
        # Robin BCs for Pathogen
        # Left boundary (x=0): v*B - d*B_x = 0
        B[0] = (4*d*B[1] - d*B[2]) / (3*d + 2*dx*v)
        
        # Right boundary (x=L): d*B_x - v*B = -c*v*B. With c=2, d*B_x = -v*B
        B[-1] = (4*d*B[-2] - d*B[-3]) / (3*d + 2*dx*v)

        # Spatial Update (Vectorized)
        dS[1:-1] = c1 * (S[2:] - 2*S[1:-1] + S[:-2]) / (dx**2) + Reaction1(S[1:-1], I[1:-1], B[1:-1], alpha_eff[1:-1], beta_eff[1:-1])
        dI[1:-1] = c2 * (I[2:] - 2*I[1:-1] + I[:-2]) / (dx**2) + Reaction2(S[1:-1], I[1:-1], B[1:-1], alpha_eff[1:-1], beta_eff[1:-1])
        dR[1:-1] = c3 * (R[2:] - 2*R[1:-1] + R[:-2]) / (dx**2) + Reaction3(I[1:-1], R[1:-1])
        dB[1:-1] = d * (B[2:] - 2*B[1:-1] + B[:-2]) / (dx**2) - v * (B[2:] - B[:-2]) / (2 * dx) + Reaction4(B[1:-1], I[1:-1], zi_eff[1:-1])
        
        # Time Step
        S += dt * dS
        I += dt * dI
        R += dt * dR
        B += dt * dB
        
    if save_history:
        S_hist.append(S.copy()); I_hist.append(I.copy())
        R_hist.append(R.copy()); B_hist.append(B.copy())
        t_hist.append(T)
        history = (np.array(S_hist), np.array(I_hist), np.array(R_hist), np.array(B_hist), np.array(t_hist))
        return S, I, R, B, history
    
    return S, I, R, B, None


# ==========================================
# 3. GRAPHING MODULE
# ==========================================
def plot_simulation(L, nx, history):
    """Takes the history output from the solver and generates 3D surface plots."""
    S_hist, I_hist, R_hist, B_hist, t_hist = history
    x = np.linspace(0, L, nx)
    X, T_mesh = np.meshgrid(x, t_hist)

    fig = plt.figure(figsize=(16, 12))
    elev_angle, azim_angle = 30, -135

    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax1.plot_surface(X, T_mesh, S_hist, cmap='turbo', edgecolor='none')
    ax1.set_title('The evolution of S'); ax1.set_xlabel('Location (x)'); ax1.set_ylabel('Time (t)')
    ax1.view_init(elev=elev_angle, azim=azim_angle)

    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    ax2.plot_surface(X, T_mesh, I_hist, cmap='turbo', edgecolor='none')
    ax2.set_title('The evolution of I'); ax2.set_xlabel('Location (x)'); ax2.set_ylabel('Time (t)')
    ax2.view_init(elev=elev_angle, azim=azim_angle)

    ax3 = fig.add_subplot(2, 2, 3, projection='3d')
    ax3.plot_surface(X, T_mesh, R_hist, cmap='turbo', edgecolor='none')
    ax3.set_title('The evolution of R'); ax3.set_xlabel('Location (x)'); ax3.set_ylabel('Time (t)')
    ax3.view_init(elev=elev_angle, azim=azim_angle)

    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    ax4.plot_surface(X, T_mesh, B_hist, cmap='turbo', edgecolor='none')
    ax4.set_title('The evolution of B'); ax4.set_xlabel('Location (x)'); ax4.set_ylabel('Time (t)')
    ax4.view_init(elev=elev_angle, azim=azim_angle)

    plt.tight_layout()
    plt.show()


# ==========================================
# 4. ACCURACY CHECK MODULE
# ==========================================
def check_spatial_accuracy(L, T, base_nx):
    """
    Implements Richardson Extrapolation to verify the spatial order of accuracy (k).
    Uses three grid resolutions: h, h/2, and h/4.
    """
    # Calculate the required dt for the FINEST grid to ensure stability across all runs
    nx_fine = 2 * (2 * base_nx - 1) - 1
    dx_fine = L / (nx_fine - 1)
    max_diffusion = max(c1, c2, c3, d)
    dt_fine = 0.4 * (dx_fine**2) / max_diffusion
    
    print(f"Running accuracy check for T={T} with fixed dt={dt_fine:.6f}...")
    
    # 1. Coarse Grid (h)
    nx_1 = base_nx
    _, I1, _, _, _ = FiniteDiffSolver(L, T, nx_1, save_history=False, fixed_dt=dt_fine)
    
    # 2. Medium Grid (h/2)
    nx_2 = 2 * nx_1 - 1
    _, I2, _, _, _ = FiniteDiffSolver(L, T, nx_2, save_history=False, fixed_dt=dt_fine)
    
    # 3. Fine Grid (h/4)
    nx_3 = 2 * nx_2 - 1
    _, I3, _, _, _ = FiniteDiffSolver(L, T, nx_3, save_history=False, fixed_dt=dt_fine)
    
    # Calculate the differences at the overlapping coarse grid points
    diff_h_h2 = I1 - I2[::2]       # w(h) - w(h/2)
    diff_h2_h4 = I2[::2] - I3[::4] # w(h/2) - w(h/4)
    
    # Calculate the L2 norm (magnitude) of these differences
    error_1 = np.linalg.norm(diff_h_h2)
    error_2 = np.linalg.norm(diff_h2_h4)
    
    # Calculate the order of accuracy k = log2( error_1 / error_2 )
    k = np.log2(error_1 / error_2)
    
    print(f"Error w(h) - w(h/2):   {error_1:.6e}")
    print(f"Error w(h/2) - w(h/4): {error_2:.6e}")
    print(f"Calculated Order of Accuracy (k): {k:.4f}")
    return k


# ==========================================
# 5. MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    L = 3 * np.pi  # Domain length from paper
    
    # --- OPTION A: Run Simulation & Graph ---
    # T = 80.0
    # nx = 100
    # print(f"Running simulation to T={T}...")
    # S, I, R, B, history = FiniteDiffSolver(L, T, nx, save_history=True)
    # plot_simulation(L, nx, history)
    
    # --- OPTION B: Run Accuracy Check ---
    T_acc = 1.0
    base_nx = 51
    check_spatial_accuracy(L, T_acc, base_nx)
