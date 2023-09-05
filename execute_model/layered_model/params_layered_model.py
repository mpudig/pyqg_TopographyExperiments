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
L = 36 * Ld      # Length of square domain

Hmax = 4000.                        # Depth of model
z = np.linspace(0, - Hmax, nz + 1)  # Cell edges 
H = z[:-1] - z[1:]                  # Height of layers (used for model)
zc = z[:-1] - H / 2                 # Cell centers



            ### Control parameters ###
    
kappa_star = 0.3
beta_star = 0.

U0 = 0.02


            ### Planetary parameters ###

g = 9.81                                  # Gravity
omega = 7.2921159e-5                      # Earth rotation angular frequency
a = 6.371e6                               # Average radius of Earth
f0 = 1e-4                                 # constant value of f0
beta = U0 / Ld ** 2 * beta_star           # Beta corresponding to the value of the nondimensional beta
lat = np.arccos(a * beta / (2 * omega))   # Latitude corresponding to that beta
# f0 = 2 * omega * np.sin(lat)            # f0 corresponding to that latitude
rek = U0 / Ld * kappa_star                # Linear friction (Ekman) drag coefficient corresponding to the nondimensional kappa



            ### Background stratification ### 
    
rho_top = 1025.                                                  # layer 1 density
rho_bot = functions.rho_bottom(f0, g, Ld, Hmax, rho_top)         # layer nz density

# Linear profile
rho = functions.linear(rho_top, rho_bot, z[:-1])

# Exponential profile
# delta = 1000.
# rho = functions.exponential(rho_top, rho_bot, z[:-1], delta)

flat_modes, flat_radii = functions.flat_bottom(f0, g, rho, H)
rough_modes, rough_radii = functions.rough_bottom(f0, g, rho, H)



            ### Background shear ###

# Project the shear onto the first baroclinic mode
BC1_mode_flat = flat_modes[:, 1] / flat_modes[0, 1]
BC1_mode_rough = rough_modes[:, 1] / rough_modes[0, 1]
U = U0 * BC1_mode_rough
V = 0 * BC1_mode_rough



            ### Random topography ###
    
K_topo = 15
h_rms = 150.
htop = functions.monoscale_random(L, nx, K_topo, h_rms)



            ### Time parameters ###

Ti = Ld / np.max(U)                  # estimate of the most unstable e-folding time scale, also nondimensionalizing factor for time [s]
dt = Ti / 400.                       # time step [s]
dt = np.ceil(dt / (60*60)) * 60*60   # time step rounded to nearest hour [s] 
tmax = 750 * Ti                      # simulation time [s]

tsnapstart = 0.                      # start time for yielding model states [s]
tsnapint = 60 * 60 * 24 * 5          # time interval at which to yield model states [s]

tavestart = 0.                       # start time to begin averaging [s]; begins halfway through the first snapshot
taveint = dt                         # time interval over which to accumulate averages [s]

tc_save = 100                               # the number of model states to store in memory to save as .nc files - should be guided by the dimensionality of the dataset         
tc_save = np.ceil(tsnapint / dt * tc_save)


            ### Initial condition ###

# The initial PV field will have energy concentrated at horizontal wavenumber K_0 and solely within vertical mode m_0, and with total energy equal to E_tot
K_0 = 9.
m_0 = 3.
q_flat_mode = flat_modes[:, int(m_0)]
q_rough_mode = rough_modes[:, int(m_0)]

# Set initial PV field
E_tot = (Ld * U.mean()) ** 2
qi = functions.set_q(K_0, q_rough_mode, L, nx, f0, g, rho, H, E_tot)


            ### Restart from previous run ###
# import glob
# import os
# restart_path = '/scratch/mp6191/pyqg_expts' + '/nx256_nz12_unforced_linearstrat_hrms150_Ktopo15_initialK24_initialm4' + '/output'
# restart_paths = sorted(glob.glob(restart_path + '/model_output_*.nc'), key = os.path.getmtime)
# restart_path = restart_paths[-1]
# qi = functions.get_q(restart_path)