import numpy as np
import xarray as xr


### Saving with snapshots or diagnostics ###
'''Note: Should use save with diagnostics if I want all the energy variables.'''

def save_with_snapshots(m, tsnapstart, tsnapint):
    '''
    Steps the model, which has already been initizlied with dt, tmax, etc., forward
    and yields (returns the model state) at an interval = tsnapint. Snapshots begin to be
    outputted at tsnapstart (units: seconds).
    
    Output:
    model_output : xarray dataset containing the model snapshots from t = 0 to t = tmax
    at a temporal resolution dt_snap = tsnapint
    '''
    
    datasets = []            # Empty list for all the model states at interval tsnapint
    m_i = m.to_dataset()     # xr dataset of initial model state
    
    for _ in m.run_with_snapshots(tsnapstart = tsnapstart, tsnapint = tsnapint):
        
        # Note: if saving averages, tsnapint should be = taveint + dt (i.e., model yields a single timestep after averages have been accumulated and saved)
        
        model_output = m.to_dataset()
        datasets.append(model_output)
        
    m_i = xr.merge([m_i, model_output]).isel(time = 0)  # Merges initial model state with final model state datasets to get diagnostic variables in initial model dataset (filled as NaNs), then removes final state dataset
    datasets.insert(0, m_i)
    
    m_ds = xr.concat(datasets,
                    dim = 'time',
                    data_vars = 'minimal')
    
    m_ds = m_ds.fillna(0.)
        
    return m_ds

def save_with_diagnostics(m, snapshots, averages, tsnapstart, tsnapint):
    '''
    Steps the model (which has already been initialized with dt, tmax, tavestart, taveint, etc.) forward
    yielding (returns the model state) at an interval tsnapint = taveint + dt.
    Then saves the averaged diagnostics (averaged over an interval taveint) chosen to be saved (the list diags).
    
    Inputs:
    m : model initialization
    snapshots : list of strings of snapshot names
    averages : list of strings of average names
    tsnapstart : the time to start yielding snapshots/returning averages at intervals
    tsnapint : the interval after which to yield snapshots/return averages
    '''
    
    datasets = []            # Empty list for all the model states at interval tsnapint
    m_i = m.to_dataset()     # xr dataset of initial model state
    m_i = m_i[snapshots]     # initial model state including only the desired diagnostics (snapshot diagnostics only)
    
    diagnostics = snapshots + averages
    
    for _ in m.run_with_snapshots(tsnapstart = tsnapstart, tsnapint = tsnapint):
                
        model_output = m.to_dataset()
        model_diags = model_output[diagnostics]
        datasets.append(model_diags)
        
    m_i = xr.merge([m_i, model_diags]).isel(time = 0)  # Merges initial model state with final model state datasets to get diagnostic variables in initial model dataset (filled as NaNs), then removes final state dataset
    datasets.insert(0, m_i)
        
    m_ds = xr.concat(datasets,
                    dim = 'time',
                    data_vars = 'minimal')
    
    return m_ds



def save_layered_model(m, snapshots, averages, tsnapstart, tsnapint, path_save, tc_save):
    '''
    Steps the layered model (which has already been initialized with dt, tmax, tavestart, taveint, etc.) forward
    yielding (returns the model state) at an interval tsnapint = taveint + dt.
    Then saves the averaged diagnostics (averaged over an interval taveint) chosen to be saved (the list diags).
    The fact that it is the layered model is important as when the time comes to save, certain attributes unique
    to the layered model have to be removed for the purposes of saving as a .nc file.
    
    Inputs:
    m : model initialization
    snapshots : list of strings of snapshot names
    averages : list of strings of average names
    tsnapstart : the time to start yielding snapshots/returning averages at intervals
    tsnapint : the interval after which to yield snapshots/return averages
    path : the path to the directory at which to save
    tc_save : the frequency with which to save (e.g., every 1000 timesteps) [units: number of model timesteps]
    '''
    
    datasets = []            # Empty list for all the model states at interval tsnapint
    m_i = m.to_dataset()     # xr dataset of initial model state
    m_i = m_i[snapshots]     # initial model state including only the desired snapshots that we will save
    
    diagnostics = snapshots + averages
    
    i = 0
    j = 0
    
    for _ in m.run_with_snapshots(tsnapstart = tsnapstart, tsnapint = tsnapint):
        model_output = m.to_dataset()
        model_diags = model_output[diagnostics]
        
        if i == 0:
            m_i = xr.merge([m_i, model_diags]).isel(time = 0)  # Merges initial model state with final model state datasets to get diagnostic variables in initial model dataset
                                                               # (filled as NaNs), then removes final state dataset
            datasets.append(m_i)
            i += 1
       
        datasets.append(model_diags)
    
        if (m.tc % tc_save) == 0:
            m_ds = xr.concat(datasets, dim = 'time', data_vars = 'minimal')   # Concatenate all datasets between timesteps j * m.tc and (j + 1) * m.tc
            del m_ds.attrs['pyqg:delta']                                      # Delete attributes that cannot be saved to a .nc file
            del m_ds.attrs['pyqg:pmodes']
            del m_ds.attrs['pyqg:radii']
            m_ds.to_netcdf(path + f'/model_output_{j}.nc')                    # Save all datasets between between timesteps j * m.tc and (j + 1) * m.tc
            del datasets                                                      # Deletes list of model states between timesteps j * m.tc and (j + 1) * m.tc
            datasets = []                                                     # Redefines empty list to be used for the next set of model states
            
            print(f'Model states between {j * m.tc + 1} and {(j + 1) * m.tc} have been saved.')
            j += 1
            
        
    m_ds = xr.concat(datasets, dim = 'time', data_vars = 'minimal')   # Concatenate all datasets between timesteps j * m.tc and the end
    del m_ds.attrs['pyqg:delta']                                      # Delete attributes that cannot be saved to a .nc file
    del m_ds.attrs['pyqg:pmodes']
    del m_ds.attrs['pyqg:radii']
    m_ds.to_netcdf(path + f'/model_output_{j + 1}.nc')                # Save all datasets between between timesteps j * m.tc and the end
            
    print(f'Model states between {j * m.tc + 1} and {m.tc} have been saved.')      
    print('Model run complete')




### Adding Qx and htop to xarray ###
'''Note: really, I want to change the native pyqg to_dataset() function to do that.'''

def to_xr_dataset(m):
    
    m_ds = m.to_dataset()
    
    # Add zonal PV gradient to xarray dataset
    Qy = m_ds.Qy
    Qx = m.Qx
    Qx = xr.DataArray(
            data = Qx,
            dims = Qy.dims,
            coords = Qy.coords,
            attrs = Qy.attrs)
    Qx.name = 'Qx'
    m_ds['Qx'] = Qx
    
    # Add bottom topography array
    
    htop = m.htop[0,:,:]
    x = m_ds.x.data
    y = m_ds.y.data
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
    m_ds['htop'] = htop
    m_ds = m_ds.fillna(0.)   # Replace nans with zeros for topography at higher levels
    
    return m_ds



### Defining topography functions ###

def band_pass(var, L, K_min, K_max):
    '''
    This function takes in a real 2D array, Fourier transforms it,
    band-pass filters between K_min and K_max,
    then inverse Fourier transforms it to return the low-pass filtered real 2D array.
    
    Inputs:
    var : real 2D array
    L : size of square domain
    K_min : normalized (i.e., integer) isotropic minimum wavenumber
    K_max : normalized (i.e., integer) isotropic maximum wavenumber 
    
    Outputs:
    var_hp : real 2D array
    '''
    
    N = var.shape[0]
    varh = np.fft.fft2(var)
    
    dk = 2 * np.pi / L
    k = dk * np.append( np.arange(0, N / 2), np.arange(- (N - 1) / 2, 0) )
    l = k
    kk, ll = np.meshgrid(k, l)
    K2 = kk ** 2 + ll ** 2
    
    K_min = K_min * (2 * np.pi / L)
    K_max = K_max * (2 * np.pi / L)
    
    # Do the band-pass on wavenumber matrix
    K2_bp_eye = np.where((K2 <= K_max ** 2) & (K2 >= K_min ** 2), K2 / K2, 0 * K2)
    K2_bp_eye[0,0] = 0
    
    # Apply filter
    varh_bp = varh * K2_bp_eye
    
    # Inverse FFT back to real space
    var_bp = np.real(np.fft.ifft2(varh_bp))
    
    return var_bp


def monoscale_random(L, N, K_min, K_max, h_rms):
    '''
    Takes in a square grid, , and minimum and maximum isotropic wavenumbers
    and returns a isotropic, homogeneous, monoscale topographic field.
    
    Inputs:
    L : side length of square box
    N : number of grid points in x and y
    K_min : normalized (i.e., integer) isotropic minimum wavenumber
    K_max : normalized (i.e., integer) isotropic maximum wavenumber
    h_rms : desired root mean square topographic height
    
    Ouputs:
    eta : topographic field (2D array shape (N, N))
    '''

    eta = band_pass(np.random.rand(N, N) * 2 - 1, L, K_min, K_max)
    
    c = h_rms / np.sqrt(np.mean(eta**2))
    
    eta = c * eta
    
    return eta


### Vertical modes ###

def flat_bottom_modes(nz, z):
    '''
    The vertical modes, \phi, the eigenfunctions of the stretching operator,
    satisfying d\phi/dz = 0 at z = 0 and d\phi/dz = 0 at z = -H
    
    NOTE: I want to also give back the corresponding eigenvalues (i.e., the deformation radii!!!) 
    '''
    
    modes = np.zeros((nz, nz))
    H = z[-1]
    
    for n in range(nz):
        modes[n,:] = np.sqrt(2) * np.cos(n * np.pi / H * z)
        
    return modes

def flat_bottom_radii(g, f0, N0, Hmax, nz):
    '''
    The flat bottom deformation radii from the eigenvalue problem. Note that this assumes linear stratification
    '''
        
    radii = np.array([N0 * Hmax / (f0 * n * np.pi) for n in range(1, nz)])  # Baroclinic deformation radii
    radii = np.insert(radii, 0, np.sqrt(g * Hmax) / f0 )                    # Barotropic deformation radii
    
    return radii


def rough_bottom_modes(nz, z):
    '''
    The vertical modes, \phi, the eigenfunctions of the stretching operator,
    satisfying d\phi/dz = 0 at z = 0 and \phi = 0 at z = -H
    
    NOTE: I want to also give back the corresponding eigenvalues (i.e., the deformation radii!!!) 
    '''
    
    modes = np.zeros((nz, nz))
    H = z[-1]
    
    for n in range(nz):
        modes[n,:] = np.sqrt(2) * np.cos( (2 * n + 1) * np.pi / (2 * H) * z)
        
    return modes

def rough_bottom_radii(f0, N0, Hmax, nz):
    '''
    The rough bottom deformation radii from the eigenvalue problem. Note that this assumes linear stratification
    '''

    radii = np.array([2 * N0 * Hmax / (f0 * (2 * n + 1) * np.pi) for n in range(1, nz)])
    
    return radii


### Setting initial PV ### 

def set_q(K_0, mode, lambda_m_0, L, nx):
    '''
    Set an initial PV distribution with energy localized in spectral space
    about K = K_0 and vertically in mode m.
    
    NOTE: ke sets the scale of the desired kinetic energy (from a scaling for U, which sets the scale of q
    
    Inputs:
    K_0 : normalized (i.e., integer) isotropic wavenumber where the energy is localized (an isotropic Gassian in normalized wavenumber space with std = 1)
    m_0 : the vertical mode number
    lambda_m_0 : the eigenvalue associated with the chosen vertical mode, i.e., the inverse deformation length squared associated with the mode
    '''
    
    # Standard deviation for |\hat \psi|^2 in isotropic wavenumber space so that std = 1 with normalized wavenumbers
    sigma = np.sqrt(2) * (2 * np.pi / L)
    
    # Define horizontal structure of PV
    dk = 2 * np.pi / L
    k = dk * np.append( np.arange(0, nx / 2), np.arange(- (nx) / 2, 0) )
    l = k
    kk, ll = np.meshgrid(k, l)
    K2 = kk ** 2 + ll ** 2
    K = np.sqrt(K2)
    K_0 = K_0 * (2 * np.pi / L)

    # Isotropic Gaussian
    psih = np.exp(-(K - K_0)**2 / (2 * sigma ** 2)) * np.exp(2 * np.pi * 1j * np.random.randn(nx, nx))
    
    # Define vertical structure of streamfunction
    psih = psih[np.newaxis, :, :] * mode[:, np.newaxis, np.newaxis]
    
    # Get qh from psih
    qh = - (K2[np.newaxis, :, :] + lambda_m_0 ** 2) * psih
    
    # Recover q from q_h, and set the scale of q given a scale for the kinetic energy
    q = np.real(np.fft.ifftn(qh, axes = (-2, -1)))
    
    ke = 1e-2
    q = q / np.max(q)                                # Normalize to have maximum unit length
    c = np.sqrt(ke) / L
    q = c * q
    
    return q


### Setting background and initial conditions ###

def linear(top, bottom, z):
    '''
    Generates a linear profile going from top to bottom with z defined positive upwards.
    '''
    
    return -(top - bottom) / z[-1] * z + top

def exponential(top, bot, z, delta):
    '''
    Generates exponential profile given a scale depth, delta,
    going from top to bottom with z defined positive upwards.
    '''
    
    return top - (top - bot) / (1 - np.exp(z[-1] / delta)) * (1 - np.exp(z/delta))

def rho_bottom(f0, g, Ld, Hmax, rho_top):
    '''
    Generates the rho_bottom which will result in a given deformation radius.
    This is, note, assuming a flat bottom deformation radius.
    '''

    const = f0 ** 2 * np.pi ** 2 * Ld ** 2 / (2 * g * Hmax)
    
    rho_bottom = (1 + const) / (1 - const) * rho_top
    
    return rho_bottom
