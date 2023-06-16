import numpy as np
import model_run_functions as functions



### Save path ###

# expt_name format: nxVALUE_nzVALUE_unforced/forced_linear/expstrat_hrmsVALUE_KtopoVALUE_K0VALUE_m0VALUE

expt_name = '/nx256_nz12_unforced_linearstrat_hrms150_Ktopo15_initialK24_initialm4'   # Modify for each experiments
path = '/scratch/mp6191/pyqg_expts' + expt_name + '/output'



### Grid ###

nx = 256         # Number of grid cells in x, y
nz = 12          # Number of grid cells in z
Ld = 50.e3       # First deformation radius
L = 16 * Ld      # Length of square domain

Hmax = 4000.                        # Depth of model
z = np.linspace(0, - Hmax, nz + 1)  # Cell edges 
H = z[:-1] - z[1:]                  # Height of layers (used for model)
zc = z[:-1] - H / 2                 # Cell centers



### Planetary parameters ###

g = 9.81                            # Gravity
f0 = 1e-4                           # Constant value of Coriolis
omega = 7.2921159e-5                # Angular frequency of Earth rotation
a = 6.371e6                         # Average radius of Earth
lat = np.arcsin(f0 / (2 * omega))   # Latitude corresponding to above f_0
beta = 2 * omega / a * np.cos(lat)  # Beta corresponding to above latitude
rek = 0.                            # Linear friction (Ekamn) drag coefficient



### Background shear, stratification, topography ###

# Zonal shear
U_top = 0.05           # layer 1 zonal velocity
U_bot = 0.             # layer nz zonal velocity

# Linear shears
U = 0 * functions.linear(U_top, U_bot, z[:-1])
V = 0 * functions.linear(U_top, U_bot, z[:-1])

# Exponential shears
# delta = 1000.                        
# U = 0 * functions.exponential(U_top, U_bot, z[:-1], delta)
# V = 0 * functions.exponential(U_top, U_bot, z[:-1], delta)

# Stratification
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

# Random bottom topography with a given rms height and with power only between two given wavenumbers
K_topo = 15
K_min = K_topo - 2
K_max = K_topo + 2
h_rms = 150.
htop = functions.monoscale_random(L, nx, K_min, K_max, h_rms)



### Time parameters ###

dt = 60 * 60 * 3                  # time step [s]
tmax = 60 * 60 * 24 * 365 * 10    # end of integration [s]

tavestart = 0.                       # start time for averaging [s]
taveint = 60 * 60 * 24 * 1           # time interval used for averaging of diagnostics [s]

tsnapstart = tavestart
tsnapint = taveint + dt



### Set initial condition for PV ###

# The initial PV field will have energy concentrated at horizontal wavenumber K_0 and solely within vertical mode m_0
K_0 = 24.
m_0 = 4.

# Flat bottom
# mode = functions.flat_bottom_modes(nz, zc)[:, int(m_0)]
# lambda_m_0 = functions.flat_bottom_radii(g, f0, N0, Hmax, nz)[int(m_0)]

# Rough bottom
mode = functions.rough_bottom_modes(nz, zc)[:, int(m_0)]
lambda_m_0 = functions.rough_bottom_radii(f0, N0, Hmax, nz)[int(m_0)]

# Set initial PV field
qi = functions.set_q(K_0, mode, lambda_m_0, L, nx)
