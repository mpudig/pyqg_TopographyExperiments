import numpy as np
import model_run_functions as functions



### Save path ###

# expt_name format: nxVALUE_nzVALUE_unforced/forced_linear/expstrat_hrmsVALUE_KtopoVALUE_K0VALUE_m0VALUE

expt_name = '/nx256_nz12_unforced_linearstrat_hrms150_Ktopo15_initialK24_initialm4'   # Modify for each experiments
path = '/scratch/mp6191/pyqg_expts' + expt_name + '/output'



            ### Grid ###

nx = 512         # Number of grid cells in x, y
nz = 12          # Number of grid cells in z
Ld = 50.e3       # First deformation radius
L = 16 * Ld      # Length of square domain

Hmax = 4000.                        # Depth of model
z = np.linspace(0, - Hmax, nz + 1)  # Cell edges 
H = z[:-1] - z[1:]                  # Height of layers (used for model)
zc = z[:-1] - H / 2                 # Cell centers



            ### Control parameters ###
    
kappa_star = 0.3
beta_star = 1.



            ### Background shear ###

# Zonal shear
U_top = 0.04           # layer 1 zonal velocity
U_bot = 0.             # layer nz zonal velocity

# Linear shears
U = 1 * functions.linear(U_top, U_bot, z[:-1])
V = 0 * functions.linear(U_top, U_bot, z[:-1])

# Exponential shears
# delta = 1000.                        
# U = 0 * functions.exponential(U_top, U_bot, z[:-1], delta)
# V = 0 * functions.exponential(U_top, U_bot, z[:-1], delta)



            ### Planetary parameters ###

g = 9.81                                  # Gravity
omega = 7.2921159e-5                      # Earth rotation angular frequency
a = 6.371e6                               # Average radius of Earth
beta = U.mean() / Ld**2 * beta_star       # Beta corresponding to the value of the nondimensional beta
lat = np.arccos(a * beta / (2 * omega))   # Latitude corresponding to that beta
# f0 = 2 * omega * np.sin(lat)              # f0 corresponding to that latitude
f0 = 1e-4                                 # constant value of f0
rek = U.mean() / Ld * kappa_star          # Linear friction (Ekman) drag coefficient corresponding to the nondimensional kappa



            ### Background stratification ### 
    
rho_top = 1025.                                                  # layer 1 density
rho_bot = functions.rho_bottom(f0, g, Ld, Hmax, rho_top)         # layer nz density

# Linear profile
rho = functions.linear(rho_top, rho_bot, z[:-1])
N = np.sqrt(- g / rho.mean() * np.gradient(rho) / np.gradient(z)[:-1])
N0 = N.mean()

# Exponential profile
# delta = 1000.
# rho = functions.exponential(rho_top, rho_bot, z[:-1], delta)
# N = np.sqrt(- g / rho.mean() * np.gradient(rho) / np.gradient(z)[:-1])


                ### Random topography ###
    
K_topo = 15
h_rms = 150.
htop = functions.monoscale_random(L, nx, K_topo, h_rms)



            ### Time parameters ###

Ti = Ld / np.abs(U_top)              # estimate of the most unstable e-folding time scale, also nondimensionalizing factor for time [s]
dt = Ti / 200.                       # time step [s]
tmax = 750 * Ti                      # simulation time [s]

tavestart = 0.                       # start time for averaging [s]
taveint = 60 * 60 * 24 * 1           # time interval used for averaging of diagnostics [s]

tsnapstart = tavestart               # start time for yielding model states
tsnapint = taveint + dt              # time interval at which to yield model states

tc_save = 100                               # the number of model states to store in memory to save as .nc files - should be guided by the dimensionality of the dataset         
tc_save = np.ceil(tsnapint / dt * tc_save)

    
# Below is what I used for the freely evolving runs:

dt = 60 * 60 * 3                     # time step [s]
tmax = 60 * 60 * 24 * 365 * 20       # end of integration [s]

tavestart = 0.                       # start time for averaging [s]
taveint = 60 * 60 * 24 * 1           # time interval used for averaging of diagnostics [s]

tsnapstart = tavestart               # start time for yielding model states
tsnapint = taveint + dt              # time interval at which to yield model states

tc_save = 100                               # the number of model states to store in memory to save as .nc files - should be guided by the dimensionality of the dataset         
tc_save = np.ceil(tsnapint / dt * tc_save)



            ### Initial condition ###

# The initial PV field will have energy concentrated at horizontal wavenumber K_0 and solely within vertical mode m_0, and with total energy equal to E_tot
K_0 = 4.
m_0 = 4.

# Flat bottom
# mode = functions.flat_bottom_modes(nz, z)[:, int(m_0)]

# Rough bottom
mode = functions.rough_bottom_modes(nz, z)[:, int(m_0)]

# Set initial PV field
qi = functions.set_q(K_0, mode, L, nx, f0, g, rho, H, E_tot)