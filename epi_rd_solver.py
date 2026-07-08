import numpy as np
import matplotlib.pyplot as plt

#Parameter Initialization
Lambda = 19.0       #population influx rate Lambda = mu * n  (Wang 2022)
                        #Lambda = .0019 * 10,000 = 19 per month
alpha = 0.000033     #base environment-to-human transmission rate (Wang 22)
beta = 0.00047      #base human to human transmission (Wang 2022)
mu = 0.0019         #natural death rate 43.5y - Wu 2024
                        #1 / (43.5 * 12) = 0.0019 per month
omega = 0.001        #disease induced death rate 
gamma = 6.0         #disease recovery rate (Range 2.9 to 14 - Wu 2024)
                        #5 days to recover: 1 / (5/30) = 6.0 per month
delta = 1.0         #removal rate of pathogen (Range 3 to 41 - Wu 2024)
                        #30 days ---> 1 / 1.0 = 1 per month
zi = 300            #base shedding rate of infected hosts (10 per day = 300 - Wang 2022)
g1 = .05            #Seasonal fluctuation amplitude of bacteria growth (Wu 2024)
g0 = 0.5           #intrinsic growth rate of the waterborne pathogen (Wu 2024)
K = 2000000.0        #carrying capacity (Wang 2022)
d = 0.1             #pathogen diffusion rate (Wang 2022)
v = 0             #pathogen advection rate 
                        #tested at v = 0, .1 , 1.2 (Wang 2022)                   


#Diffusion rates (Wang 2022)
c1 = 1.0  #succeptible 
c2 = 0.5  #infected
c3 = 0.8  #recovered

def Reaction1(S, I, B, alpha_t, beta_t):
    return Lambda - alpha_t*S*B - beta_t*S*I - mu*S

def Reaction2(S, I, B, alpha_t, beta_t):
    return alpha_t*S*B + beta_t*S*I - (mu + omega + gamma)*I

def Reaction3(I, R):
    return gamma*I - mu*R

def Reaction4(B, I, zi_t, g_t):
    return g_t*B * (1 - B/K) + zi_t*I - delta*B


#Numerical Method Solver FTCS
def FiniteDiffSolver(L, T, nx):  
    """
    Solves the 4-equation SIR-B model using FTCS finite difference.
    Stores the entire history in a 2D array of size (nt, nx).
    """
    dx = L / (nx - 1)

    #Calculate stable time step 
    max_diffusion = max(c1, c2, c3, d)
    dt_max = 0.4 * (dx**2) / max_diffusion
    nt = int(T / dt_max) + 1

    
    #time array
    dt = T / nt
    t_array = np.linspace(0, T, nt)

    #Initialize 2D arrays
    S = np.zeros((nt, nx))
    I = np.zeros((nt, nx))
    R = np.zeros((nt, nx))
    B = np.zeros((nt, nx))

    #Spatial grid for initial conditions
    x = np.linspace(0, L, nx)
    H_x = 0.5 + 0.25 * np.cos(2 * x)

    #Initial Conditions (Wang 2022) pg 16
    S[0, :] = (Lambda / mu) - 500 * np.cos(2 * x) 
    I[0, :] = 1 - np.cos(2 * x)
    R[0, :] = np.zeros(nx)
    B[0, :] = 0.5 - 0.3 * np.cos(2 * x)
    
    #Temporal Loop 
    for n in range(nt - 1):
        #Seasonality Function
        T_t = 0.5 + 0.25 * np.sin(2 * np.pi * t_array[n] / 12.0) #Wang 2022 pg 15
        
        # m_I = I[n, :] / (I[n, :] + M) # Wang 2022 pg 16
        # alpha_t = alpha * (1 - b_factor * m_I) * H_x * T_t 
        # beta_t = beta * (1 - b_factor * m_I) * H_x * T_t #Wu 2024
        # zi_t = zi * (1 - b_factor * m_I) * H_x * T_t #Wang 2022
        # g_t = g0 + 0.5 * np.sin(np.pi * t_array[n] / 6.0) #Wang 2022 

        #Effective rates
        alpha_t = alpha * H_x * T_t
        beta_t = beta * H_x * T_t #Wu 2024 pg 30
        zi_t = zi * H_x * T_t #Wang 2022
        g_t = g0 + g1 * np.sin(np.pi * t_array[n] / 6.0) #Wang 2022 / Wu 2024 (pg 15)

        #Spatial Update for interior nodes (1 to nx-2)
        dS = c1 * (S[n, 2:] - 2*S[n, 1:-1] + S[n, :-2]) / (dx**2) + Reaction1(S[n, 1:-1], I[n, 1:-1], B[n, 1:-1], alpha_t[1:-1], beta_t[1:-1])
        dI = c2 * (I[n, 2:] - 2*I[n, 1:-1] + I[n, :-2]) / (dx**2) + Reaction2(S[n, 1:-1], I[n, 1:-1], B[n, 1:-1], alpha_t[1:-1], beta_t[1:-1])
        dR = c3 * (R[n, 2:] - 2*R[n, 1:-1] + R[n, :-2]) / (dx**2) + Reaction3(I[n, 1:-1], R[n, 1:-1])
        dB = d * (B[n, 2:] - 2*B[n, 1:-1] + B[n, :-2]) / (dx**2) - v * (B[n, 2:] - B[n, :-2]) / (2 * dx) + Reaction4(B[n, 1:-1], I[n, 1:-1], zi_t[1:-1], g_t)
        
        #Time Step Update
        S[n+1, 1:-1] = S[n, 1:-1] + dt * dS
        I[n+1, 1:-1] = I[n, 1:-1] + dt * dI
        R[n+1, 1:-1] = R[n, 1:-1] + dt * dR
        B[n+1, 1:-1] = B[n, 1:-1] + dt * dB

        #Apply Boundary Conditions to n+1 row
        S[n+1, 0] = (4*S[n+1, 1] - S[n+1, 2]) / 3.0
        S[n+1, -1] = (4*S[n+1, -2] - S[n+1, -3]) / 3.0
        
        I[n+1, 0] = (4*I[n+1, 1] - I[n+1, 2]) / 3.0
        I[n+1, -1] = (4*I[n+1, -2] - I[n+1, -3]) / 3.0
        
        R[n+1, 0] = (4*R[n+1, 1] - R[n+1, 2]) / 3.0
        R[n+1, -1] = (4*R[n+1, -2] - R[n+1, -3]) / 3.0
        
        #Apply Pathogen Boundary Conditions
        B[n+1, 0] = (4*d*B[n+1, 1] - d*B[n+1, 2]) / (3*d + 2*dx*v)
        B[n+1, -1] = (4*d*B[n+1, -2] - d*B[n+1, -3]) / (3*d + 2*dx*v)
        
    
    return S, I, R, B, t_array



def graphResults(L, nx, S, I, R, B, t_array):
    """
    Plots the 2D arrays from the numerical method as 3D graphs
    """
    x = np.linspace(0, L, nx)
    X, T_mesh = np.meshgrid(x, t_array)

    fig = plt.figure(figsize=(10, 8))
    elev_angle, azim_angle = 30, -135

    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax1.plot_surface(X, T_mesh, S, cmap='turbo', edgecolor='none')
    ax1.set_title('Susceptible Individuals (S)'); ax1.set_xlabel('Location (x)'); ax1.set_ylabel('Time (t)')
    ax1.view_init(elev=elev_angle, azim=azim_angle)

    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    ax2.plot_surface(X, T_mesh, I, cmap='turbo', edgecolor='none')
    ax2.set_title('Infected Individuals (I)'); ax2.set_xlabel('Location (x)'); ax2.set_ylabel('Time (t)')
    ax2.view_init(elev=elev_angle, azim=azim_angle)

    ax3 = fig.add_subplot(2, 2, 3, projection='3d')
    ax3.plot_surface(X, T_mesh, R, cmap='turbo', edgecolor='none')
    ax3.set_title('Recovered Individuals (R)'); ax3.set_xlabel('Location (x)'); ax3.set_ylabel('Time (t)')
    ax3.view_init(elev=elev_angle, azim=azim_angle)

    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    ax4.plot_surface(X, T_mesh, B, cmap='turbo', edgecolor='none')
    ax4.set_title('Pathogen (B)'); ax4.set_xlabel('Location (x)'); ax4.set_ylabel('Time (t)')
    ax4.view_init(elev=elev_angle, azim=azim_angle)

    fig.supxlabel(f'v = {v}', fontsize=16)
    plt.tight_layout()
    plt.show()


def spatial_accuracy_check(L, T, nx):
    """
    This function checks the spatial order of accuracy for the method.
    """
 
    #Calculate w(h)
    nx_1 = nx
    _, I1, _, _, _ = FiniteDiffSolver(L, T, nx_1)
    
    #Calculate w(h/2)
    nx_2 = 2 * nx_1 - 1
    _, I2, _, _, _ = FiniteDiffSolver(L, T, nx_2)
    
    #Calculate w(h/4)
    nx_3 = 2 * nx_2 - 1
    _, I3, _, _, _ = FiniteDiffSolver(L, T, nx_3)
    
    #Calculate the differences 
    diff_h_h2 = I1[-1, :] - I2[-1, ::2]       # w(h) - w(h/2)
    diff_h2_h4 = I2[-1, ::2] - I3[-1, ::4]    # w(h/2) - w(h/4)
    
    #Calculate the L2 norm of differences 
    error_1 = np.linalg.norm(diff_h_h2)
    error_2 = np.linalg.norm(diff_h2_h4)
    
    #Calculate the order of accuracy 
    k = np.log2(error_1 / error_2)
    
    print(f"Error w(h) - w(h/2):   {error_1:.6e}")
    print(f"Error w(h/2) - w(h/4): {error_2:.6e}")
    print(f"Calculated Order of Accuracy (k): {k:.4f}")

    #calculate dx values
    dx_1 = L / (nx_1 - 1)
    dx_2 = L / (nx_2 - 1)
    
    dx_values = np.array([dx_1, dx_2])
    errors = np.array([error_1, error_2])
    
    #Create reference line of slope 2
    reference_line = errors[0] * (dx_values / dx_1)**2 

    #log-log convergence plot
    plt.figure(figsize=(8, 6))
    plt.loglog(dx_values, errors, 'bo-', linewidth=2, markersize=8, label='Numerical Error $||w_h - w_{h/2}||$')
    plt.loglog(dx_values, reference_line, 'k--', linewidth=2, label=r'Reference Slope = 2 ($\mathcal{O}(\Delta x^2)$)')
    
    plt.xlabel(r'Spatial Step Size ($\Delta x$)')
    plt.ylabel('L2 Norm of Error')
    plt.title('Spatial Grid Convergence')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.show()

    return k

if __name__ == "__main__":
    L = 3 * np.pi
    T = 12
    nx = 120
    
    S, I, R, B, t_array = FiniteDiffSolver(L, T, nx)
    
    #graphResults(L, nx, S, I, R, B, t_array)

    spatial_accuracy_check(L, T, nx)

