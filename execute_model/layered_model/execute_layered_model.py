import params_layered_model as params
import diags_layered_model as diags
import model_run_functions as functions
import pyqg
import numpy as np
import xarray as xr

            ### Grid ###

nx = params.nx
nz = params.nz
L = params.L
H = params.H

            ### Planetary parameters ###

f0 = params.f0
beta = params.beta
rek = params.rek

            ### Background shear, stratification, topography ###

U = params.U
V = params.V
rho = params.rho
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

m = pyqg.LayeredModel(nx = nx, nz = nz, L = L, H = H,
                      U = U, V = V, rho = rho, htop = htop,
                      f = f0, beta = beta, rek = rek,
                      dt = dt, tmax = tmax, twrite = 10, tavestart = tavestart, taveint = taveint, ntd = ntd)

            ### Set initial condition for PV ###

m.set_q(params.qi)



            ### Run model with diagnostics at given frequency, and save at given path ###

snapshots = diags.snapshots
averages = diags.averages
path = params.path
tc_save = params.tc_save

# Save topography field
functions.save_htop(m, htop, path)

# Run and save model
functions.save_layered_model(m, snapshots, averages, tsnapstart, tsnapint, path, tc_save)
