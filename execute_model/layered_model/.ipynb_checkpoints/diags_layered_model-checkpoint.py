### Available diagnostics (snapshots and averages) for layered model ###

# Commment out those not wanted #

snapshots = [
    'q',                   # potential vorticity in real space
    'p',                   # streamfunction in real space
    'u',                   # zonal velocity anomaly
    'v',                   # meridional velocity anomaly
#     'ufull',               # zonal full velocities in real space
#     'vfull',               # meridional full velocities in real space
#     'Ubg',                 # background zonal velocity
#     'Qy',                  # background potential vorticity gradient
#     'qh',                  # potential vorticity in spectral space
#     'uh',                  # zonal velocity anomaly in spectral space
#     'vh',                  # meridional velocity anomaly in spectral space
#     'ph',                  # streamfunction in spectral space
#     'dqdt',                # previous partial derivative of potential vorticity wrt. time in real space
#     'dqhdt',               # previous partial derivative of potential vorticity wrt. time in spectral space
]

# Reorder the below for ease of readability, lots of them I won't use:
averages = [
#     'EKE',                 # mean eddy kinetic energy
#     'KEspec',              # kinetic energy spectrum
#     'APEspec',             # available potential energy spectrum 
#     'Ensspec',             # enstrophy spectrum
#     'KEspec_modal',        # modal kinetic energy spectra (NOTE: These are currently _flat_ modes)
#     'PEspec_modal',        # modal potential energy spectra (NOTE: These are currently _flat_ modes)
#     'entspec',             # barotropic enstrophy spectrum
#     'APEgenspec',          # the spectrum of the rate of generation of available potential energy
#     'Dissspec',            # spectral contribution of filter dissipation to total energy
#     'KEfrictionspec',      # total energy dissipation spectrum by bottom drag
#     'KEflux_div',          # spectral divergence of flux of kinetic energy
#     'APEflux_div',         # spectral divergence of flux of available potential energy
#     'EKEdiss',             # total energy dissipation by bottom drag
#     'ENSDissspec',         # spectral contribution of filter dissipation to barotropic enstrophy
#     'ENSflux',             # barotropic enstrophy flux
#     'ENSfrictionspec',     # the spectrum of the rate of dissipation of barotropic enstrophy due to bottom friction
#     'ENSgenspec',          # the spectrum of the rate of generation of barotropic enstrophy
#     'paramspec',           # spectral contribution of subgrid parameterization to energy (if present)
#     'ENSparamspec',        # spectral contribution of subgrid parameterization to enstrophy (if present)
#     'paramspec_KEflux',    # total additional KE flux due to subgrid parameterization   
#     'paramspec_APEflux',   # total additional APE flux due to subgrid parameterization
]