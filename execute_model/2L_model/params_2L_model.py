import numpy as np
import model_run_functions as functions



            ### Save path ###

expt_name = '/<expt_name>'                                   # Modify for each experiments
path = '/scratch/mp6191/pyqg_expts' + expt_name + '/output'  # Path for saving model output



            ### Grid ###

nx = 512                            # Number of grid cells in x, y
Ld = 15.e3                          # deformation radius
kd = 1 / Ld                         # deformation wavenumber
ld = 2 * np.pi * Ld                 # deformation wavelength
L = 25 * ld                         # length of domain

H = 4000.                           # total depth
delta = 1.                          # layer ratio
H1 = delta / (1 + delta) * H        # height of first layer
H2 = 1 / (1 + delta) * H            # height of second layer   
Hi = np.array([H1, H2])             # array of layers
htop = np.zeros((nx, nx))           # topography (need to fix model blow-up when topography not defined!)


            ### Control parameters ###
    
kappa_star = 0.3                    # kappa* = kappa * ld / U = 1/(2*pi) * kappa * Ld / U
beta_star = 0.6                     # beta* = beta * ld^2 / U = 1/(2*pi) * beta * Ld^2 / U

U = 0.01                            # mean value of background velocity


            ### Planetary parameters ###

g = 9.81                                  # Gravity
omega = 7.2921159e-5                      # Earth rotation angular frequency
a = 6.371e6                               # Average radius of Earth
# f0 = 1e-4                                # constant value of f0
beta = U / Ld ** 2 * beta_star            # Beta corresponding to the value of the nondimensional beta
lat = np.arccos(a * beta / (2 * omega))   # Latitude corresponding to that beta
f0 = 2 * omega * np.sin(lat)              # f0 corresponding to that latitude
rek = U / Ld * kappa_star                 # linear friction (Ekman) drag coefficient corresponding to the nondimensional kappa



            ### Background shear ###
    
U1 = 2 * U                          # upper layer imposed flow
U2 = 0.                             # lower layer imposed flow


            ### Time parameters and threading ###

Ti = Ld / U                  # estimate of the most unstable e-folding time scale, also nondimensionalizing factor for time [s]
dt = 3600.                           # time step [s]
tmax = 750 * Ti                      # simulation time [s]

tsnapstart = 0.                      # start time for yielding model states [s]
tsnapint = 60 * 60 * 24              # time interval at which to yield model states [s]

tavestart = 0.                       # start time to begin averaging [s]; begins halfway through the first snapshot
taveint = dt                         # time interval over which to accumulate averages [s]

tc_save = 50                         # the number of model states to store in memory to save as .nc files - should be guided by the dimensionality of the dataset         
tc_save = np.ceil(tsnapint / dt * tc_save)

ntd = 16


            ### Initial condition ###

# The initial PV field will have energy concentrated at horizontal wavenumber K_0 and solely within vertical mode m_0, and with total energy equal to E_tot
K_0 = L / (4 * Ld)
E_tot = (Ld * U) ** 2

# qi = functions.set_q(K_0, L, nx, f0, g, Hi, E_tot)


            ### Restart from previous run ###
# import glob
# import os
# restart_path = '/scratch/mp6191/pyqg_expts' + '/nx256_nz12_unforced_linearstrat_hrms150_Ktopo15_initialK24_initialm4' + '/output'
# restart_paths = sorted(glob.glob(restart_path + '/model_output_*.nc'), key = os.path.getmtime)
# restart_path = restart_paths[-1]
# qi = functions.get_q(restart_path)