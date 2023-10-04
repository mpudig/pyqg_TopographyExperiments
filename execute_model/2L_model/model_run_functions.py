import numpy as np
import xarray as xr

            ### Saving with snapshots or diagnostics ###

def save_2L_model(m, snapshots, averages, tsnapstart, tsnapint, path, tc_save):
    '''
    This function steps the 2L model (which has already been initialized with dt, tmax, tavestart, taveint, etc.)
    forward yielding (returns the model state) at an interval tsnapint. The model state snapshot should include the desired
    instantenous snapshot variables (e.g., q) as well as the averaged diagnostics that have been accumulating with a
    frequency taveint.
    
    
    Inputs:
    m : model initialization
    snapshots : list of strings of snapshot names
    averages : list of strings of average names
    tsnapstart : the time to start yielding snapshots/returning averages at intervals
    tsnapint : the interval after which to yield snapshots/return averages
    path : the path to the directory at which to save
    tc_save : the frequency with which to save (e.g., every 1000 timesteps) [units: number of model timesteps]
    '''
    
    datasets = []                       # Empty list for all the model states at interval tsnapint
    m_i = (m.to_dataset())[snapshots]   # xr dataset of initial model state

    # Set all snapshots and averages we want to save
    diagnostics = snapshots + averages

    # Set other averages as inactive so that pyqg doesn't calculate them as well, to reduce computational burden
    averages_all = list(m.diagnostics.keys())
    averages_inactive = [average for average in averages_all if average not in averages]

    for d in list(averages_inactive):
        m.diagnostics[d]['active'] = False

    
    i = 0
    j = 0
    m_tc_init = 0

    # Run the model forward
    for _ in m.run_with_snapshots(tsnapstart = tsnapstart, tsnapint = tsnapint):
        
        # Make xarray dataset including desired snapshots and averages 
        model_output = (m.to_dataset())[diagnostics]
        
        # Properly add the initial condition model state to the datasets list, with average diagnostics set to np.nan   
        if i == 0:
            m_i = xr.merge([m_i, model_output]).isel(time = 0)
            datasets.append(m_i)
            i += 1
        
        # Append the xarray dataset model state to the list of datasets
        datasets.append(model_output)
        
        # Reset all average diagnostics so that the new set of averages begins and is available at next snapshot
        m._initialize_diagnostics('all')
        for d in list(averages_inactive):
            m.diagnostics[d]['active'] = False

        # Save as .nc file    
        if (m.tc % tc_save) == 0:
            m_ds = xr.concat(datasets, dim = 'time', data_vars = 'minimal')   # Concatenate all datasets between given timesteps
            
            # Delete annoying attributes which stuff up the saving
            del m_ds.attrs['pyqg:delta']                                      # Delete attributes that cannot be saved to a .nc file
            del m_ds.attrs['pyqg:pmodes']
            del m_ds.attrs['pyqg:radii']
            
            # Save to path and then delete the datasets from memory
            m_ds.to_netcdf(path + f'/model_output_{j}.nc')                    # Save all datasets between between given timesteps
            del datasets                                                      # Deletes list of model states between given timesteps
            datasets = []                                                     # Redefines empty list to be used for the next set of model states
            
            # Tell me what happened
            print(f'Model states between {m_tc_init} and {m.tc} have been saved.')
            m_tc_init = m.tc
            j += 1
            
        
    m_ds = xr.concat(datasets, dim = 'time', data_vars = 'minimal')   # Concatenate all datasets between given timesteps and the end

    m_ds.to_netcdf(path + f'/model_output_{j}.nc')                # Save all datasets between given timesteps and end
            
    print(f'Model states between {m_tc_init + 1} and {m.tc} have been saved.')      
    print('Model run complete')

    
    
            ### Saving topography field ###

def save_htop(m, htop, path):
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



            ### Defining topography functions ###

def monoscale_random(L, N, K_0, h_rms):
    '''
    Takes in a square grid, and minimum and maximum isotropic wavenumbers
    and returns a isotropic, homogeneous, monoscale topographic field.
    
    Inputs:
    L : side length of square box
    N : number of grid points in x and y
    K_0 : central normalized (i.e., integer) isotropic minimum wavenumber
    h_rms : desired root mean square topographic height
    
    Ouputs:
    eta : topographic field (2D array shape (N, N))
    '''
    
    # Standard deviation for |\hat \psi|^2 in isotropic wavenumber space so that std = 1 with normalized wavenumbers
    sigma = np.sqrt(2) * (2 * np.pi / L)
    
    # Define horizontal structure of PV
    dk = 2 * np.pi / L
    k = dk * np.append( np.arange(0, N / 2), np.arange(- N / 2, 0) )
    l = k
    kk, ll = np.meshgrid(k, l)
    K2 = kk ** 2 + ll ** 2
    K = np.sqrt(K2)
    K_0 = K_0 * (2 * np.pi / L)
    
    # Isotropic Gaussian in wavenumber space
    etah = np.exp(-(K - K_0) ** 2 / (2 * sigma ** 2)) * np.exp(2 * np.pi * 1j * np.random.randn(N, N))

    # Recover eta from eta_h
    eta = np.real(np.fft.ifftn(etah, axes = (-2, -1)))
    
    c = h_rms / np.sqrt(np.mean(eta**2))
    eta = c * eta
    
    return eta


            ### Vertical modes and radii ###
    
def flat_bottom(f0, g, rho, H):
    '''
    Calculates the flat bottom vertical modes and radii (i.e., those associated with zero flux BCs at the surface and at the bottom).
    '''
    
    ### Generate flat bottom stretching matrix
    nz = H.size
    S = np.zeros((nz, nz))
    
    f2 = f0 ** 2
    gpi = g * (rho[1:] - rho[:-1]) / rho[:-1]
    
    for i in range(nz):
        
        if i == 0:
            S[i, i]   = - f2 / H[i] / gpi[i]
            S[i, i + 1] =  f2 / H[i] / gpi[i]

        elif i == nz - 1:
            S[i, i]   = - f2 / H[i] / gpi[i - 1]
            S[i, i - 1] =  f2 / H[i] / gpi[i - 1]

        else:
            S[i, i - 1] = f2 / H[i] / gpi[i - 1]
            S[i, i] = - (f2 / H[i] / gpi[i] + f2 / H[i] / gpi[i - 1])
            S[i, i + 1] = f2 / H[i] / gpi[i]
            
            
    ### Calculate eigenvectors and eigenvalues of matrix
        
    evals, evecs = np.linalg.eig(-S)
    asort = evals.argsort()
    
    # deformation radii
    evals = evals[asort]
    radii = 1. / np.sqrt(evals)
    radii[0] = np.sqrt(g * H.sum()) / np.abs(f0)  # barotropic def. radius
    
    # eigenstructure
    modes = evecs[:, asort]
    
    # normalize to have unit L2-norm
    Ai = (H.sum() / (H[:, np.newaxis] * (modes**2)).sum(axis=0)) ** 0.5
    modes = Ai[np.newaxis, :] * modes
    
    return modes, radii


def rough_bottom(f0, g, rho, H):
    '''
    Calculates the rough bottom vertical modes and radii (i.e., those associated with zero flux BC at the surface and zero Dirichlet at the bottom).
    One has to extrapolate the rho array to have a value in the layer below that explicitly modelled in order to have second order accuracy at the
    bottom boundary. Note that the way I currently have this set up assumes that the rho profile is linear. I should modify this for when I will have
    exponential, or other more general, stratifications.
    '''
    
    ### Generate rough bottom stretching matrix
    nz = H.size
    S = np.zeros((nz, nz))
    
    f2 = f0 ** 2
    # Extrapolate density to have a value in the layer below that modelled. Note that this currently assumes a linear density profile.
    drho = np.diff(rho)[0]
    rho = np.append(rho, rho[-1] + drho)
    gpi = g * (rho[1:] - rho[:-1]) / rho[:-1]
    
    for i in range(nz):
        
        if i == 0:
            S[i, i]   = - f2 / H[i] / gpi[i]
            S[i, i + 1] =  f2 / H[i] / gpi[i]

        elif i == nz - 1:
            S[i, i]   = - f2 / H[i] / gpi[i - 1]
            S[i, i - 1] =  f2 / H[i] / gpi[i - 1]

        else:
            S[i, i - 1] = f2 / H[i] / gpi[i - 1]
            S[i, i] = - (f2 / H[i] / gpi[i] + f2 / H[i] / gpi[i - 1])
            S[i, i + 1] = f2 / H[i] / gpi[i]
            
    # Rough bottom stretching matrix simply adds a term to the S_nn entry of the flat bottom stretching matrix
    S[nz - 1, nz - 1] -= 2 * f2 / (H[-1] * gpi[-1])
            
            
    ### Calculate eigenvectors and eigenvalues of matrix
        
    evals, evecs = np.linalg.eig(-S)
    asort = evals.argsort()
    
    # deformation radii
    evals = evals[asort]
    radii = 1. / np.sqrt(evals)
    
    # eigenstructure
    modes = evecs[:, asort]
    
    # normalize to have unit L2-norm
    Ai = (H.sum() / (H[:, np.newaxis] * (modes**2)).sum(axis=0)) ** 0.5
    modes = Ai[np.newaxis, :] * modes
    
    return modes, radii
    

            ### Setting initial PV ### 

def set_q(K_0, eigenvector, L, nx, f0, g, rho, H, E_tot):
    '''
    Set an initial PV distribution with energy localized in spectral space
    about K = K_0 and vertically in mode m. The PV field is also rescaled so that
    the total energy is (approximately, it depends on the vertical discretization unfortunately)
    given by E_tot.
    
    Inputs:
    K_0 : normalized (i.e., integer) isotropic wavenumber where the energy is localized (an isotropic Gassian in normalized wavenumber space with std = 1)
    eigenvector : the vertical mode that sets the vertical structure of q
    L : the length of the horizontal domain
    nx : the number of grid cells in the horizontal domain
    f0 : constant value of Coriolis
    g : gravity
    rho : array of imposed background densities in each layer
    H : array of average depths in each layer
    E_tot : the prescribed total energy for the initial condition
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

    # \hat \psi is isotropic Gaussian, we also want it to correspond to a real field so force it to do so
    psih = np.exp(-(K - K_0)**2 / (2 * sigma ** 2)) * np.exp(2 * np.pi * 1j * np.random.randn(nx, nx))
    psi = np.fft.irfftn(psih[:, :int(nx / 2) + 1], axes = (-2,-1))
    psih = np.fft.fftn(psi, axes = (-2, -1))
    
    # Define vertical structure of streamfunction to have baroclinic structure
    psih = psih[np.newaxis, :, :] * np.array([1, -1])[:, np.newaxis, np.newaxis]
    
    # Calculate total energy and scaling factor so that energy of nondimensional system is unity
    f2 = f0 ** 2
    gpi = g * (rho[1:] - rho[:-1]) / rho[:-1]   # Reduced gravities
    M = nx * nx                                 # Spectral normalization
    
    KE = L ** 2 * 1 / (2 * H.sum()) * (H[:, np.newaxis, np.newaxis] * K2[np.newaxis, :, :] * np.abs(psih / M)**2).sum()
    APE = L ** 2 * 1 / (2 * H.sum()) * (f2 / gpi[:, np.newaxis, np.newaxis] * np.abs(psih[:-1, :, :] / M - psih[1:, :, :] / M) ** 2).sum()
    E = KE + APE
    c = np.sqrt(E_tot / E)
    psih = c * psih

    
    # Calculate matrix to get qh from psih
    nz = H.size
    
    S = np.zeros((nz, nz))
    for i in range(nz):
        
        if i == 0:
            S[i, i]   = - f2 / H[i] / gpi[i]
            S[i, i + 1] =  f2 / H[i] / gpi[i]

        elif i == nz - 1:
            S[i, i]   = - f2 / H[i] / gpi[i - 1]
            S[i, i - 1] =  f2 / H[i] / gpi[i - 1]

        else:
            S[i, i - 1] = f2 / H[i] / gpi[i - 1]
            S[i, i] = - (f2 / H[i] / gpi[i] + f2 / H[i] / gpi[i - 1])
            S[i, i + 1] = f2 / H[i] / gpi[i]
            
    I = np.eye(nz)[:, :, np.newaxis, np.newaxis]
    M = S[:, :, np.newaxis, np.newaxis] - I * K2
    
    # Get qh from psih
    qh = np.zeros_like(psih)
    for m1 in range(nz):
        for m2 in range(nz):
            for j in range(nx):
                for i in range(nx):
                    qh[m2, j, i] = qh[m2, j, i] + M[m2, m1, j, i] * psih[m1, j, i]
        
    # Recover q    
    q = np.real(np.fft.ifftn(qh, axes = (-2, -1)))
    
    return q


def get_q(path):
    '''                                                                                                                                    
    Picks up the last snapshot of q from a previous run to be used as a restart field for a new run.
                                                                                                                                           
    path : the complete path to the .nc file which contains q.                                                                             
    '''

    arr = xr.open_dataset(path)
    q = arr.q.isel(time = -1).load().values

    return q


            ### Setting background and initial conditions ###
    
def rho_bottom(f0, g, Ld, Hmax, rho_top):
    '''
    Generates the rho_bottom which will result in a given deformation radius.
    This is, note, assuming a flat bottom deformation radius.
    '''

    const = f0 ** 2 * np.pi ** 2 * Ld ** 2 / (2 * g * Hmax)
    
    rho_bottom = (1 + const) / (1 - const) * rho_top
    
    return rho_bottom

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