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

### Time parameters ###

dt = params.dt
tmax = params.tmax

tavestart = params.tavestart
taveint = params.taveint

tsnapstart = params.tsnapstart
tsnapint = params.tsnapint

### Initialize model instance ###

m = pyqg.LayeredModel(nx = nx, nz = nz, L = L, H = H,
                      U = U, V = V, rho = rho, htop = htop,
                      f = f0, beta = beta, rek = rek,
                      dt = dt, tmax = tmax, tavestart = tavestart, taveint = taveint)

### Set initial condition for PV ###

m.set_q(params.qi)



### Run model with diagnostics at given frequency, and save at given path ###

snapshots = diags.snapshots
averages = diags.averages
path = params.path

# Save topography field 
x = m.x[0, :]
y = m.y[:, 0]
htop = xr.DataArray(
            data = htop,
            dims = ['x', 'y'],
            coords = {
                'x': ('x', x),
                'y': ('y', y)},
            attrs = dict(
                units = 'm',
                long_name = 'height of bottom topography field'))
htop.name = 'htop'
htop_path = path + '/htop.nc'
htop.to_netcdf(htop_path)

# Run and save model
m_ds = functions.save_with_diagnostics(m, snapshots, averages, tsnapstart, tsnapint)

# Removes attributes for saving purposes
del m_ds.attrs['pyqg:delta']
del m_ds.attrs['pyqg:pmodes']
del m_ds.attrs['pyqg:radii']

# Save model output to path
m_ds.to_netcdf(path + '/model_output.nc')
print('Run and saving complete')