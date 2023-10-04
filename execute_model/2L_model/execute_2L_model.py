import params_2L_model as params
import diags_2L_model as diags
import model_run_functions as functions
import pyqg
import numpy as np
import xarray as xr

            ### Grid ###

nx = params.nx
L = params.L
H1 = params.H1

            ### Planetary parameters ###

f0 = params.f0
beta = params.beta
rek = params.rek

            ### Background shear, stratification, topography ###

U1 = params.U1
U2 = params.U2
htop = params.htop

            ### Time parameters and threading ###

dt = params.dt
tmax = params.tmax

tavestart = params.tavestart
taveint = params.taveint

tsnapstart = params.tsnapstart
tsnapint = params.tsnapint

ntd = params.ntd

            ### Initialize model instance ###

m = pyqg.QGModel(nx = nx, L = L, H1 = H1, U1 = U1, U2 = U2, f = f0, beta = beta, rek = rek, htop = htop,
                       dt = dt, tmax = tmax, twrite = 100)

# m = pyqg.QGModel_rough(nx = nx, L = L, H1 = H1, U1 = U1, U2 = U2, f = f0, beta = beta, rek = rek, htop = htop,
#                        dt = dt, tmax = tmax, twrite = 100)

            ### Set initial condition for PV ###

# m.set_q(params.qi)



            ### Run model with diagnostics at given frequency, and save at given path ###

snapshots = diags.snapshots
averages = diags.averages
path = params.path
tc_save = params.tc_save

# Run and save model
functions.save_2L_model(m, snapshots, averages, tsnapstart, tsnapint, path, tc_save)