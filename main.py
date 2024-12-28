from eady import eady_model
import numpy as np 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from functions import inverse_transform, wavenumbers, transform
from matplotlib.animation import FuncAnimation
from scipy.stats import linregress

# USUAL PARAMETERS
Lx = 8e6;                   Ly = Lx;                 	Lz = 1e4
Nx = 2**6;                  Ny = 2**4;                  Nz = 50
tmax = 3600*24*30;          dt=3900;                    Umax = 10

N = 0.01
x = np.linspace(0, Nx, Nx)
y = np.linspace(0, Ny, Ny)
X, Y = np.meshgrid(x, y)
z = np.linspace(0,Lz,Nz)

def eady_analysis(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45, f_cte=True, perturbation=2, linear_model=False, rho_cte=True, custom_profile=False, U0=None, dU0dz=None, d2U0dz2=None):
    """
    Perform Eady model analysis to compute the growth rate and other related parameters.
    
    Outputs:
        - growth_rate (float): Growth rate of the instability.
        - times (numpy.ndarray): Array of time values.
        - maxV_values (numpy.ndarray): Array of maximum velocity values.
        - Q (numpy.ndarray): Array of Q values.
        - v (numpy.ndarray): Array of v values.
        - U (numpy.ndarray): Array of U values.
        - PSI (numpy.ndarray): Array of streamfunction values.
        - zvorticity (numpy.ndarray): Array of vorticity values.
    """

    time, maxV_values, Q, v, U, PSI, zvorticity = eady_model(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude, f_cte,perturbation, linear_model, rho_cte, custom_profile=custom_profile, U0=U0, dU0dz=dU0dz, d2U0dz2=d2U0dz2)

    times = np.array(time)
    Omega = 7.2921e-5
    f0 = 2*Omega*np.sin(np.deg2rad(latitude))
    Ld=N*Lz/f0
    ss = times/Ld*Umax
    ss = ss[int(0.8 * len(ss)):]
    log_maxV = np.log(maxV_values)
    log_maxV = log_maxV[int(0.8 * len(log_maxV)):]
    growth_rate = 0
    if tmax > 0:  
        slope, intercept, r_value, p_value, std_err = linregress(ss, log_maxV)
        growth_rate = slope  
        
    return growth_rate, times, maxV_values, Q, v, U, PSI, zvorticity

########################################################################################

def plot_growth_rate_setups(setups, Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45):
    """
    Plots the growth rates of different setups using the Eady model analysis.

    Parameters:
    setups (list of dict): A list of dictionaries, each containing the parameters for a specific setup.
    Lx (float): Length of the domain in the x-direction.
    Ly (float): Length of the domain in the y-direction.
    Lz (float): Length of the domain in the z-direction.
    Nx (int): Number of grid points in the x-direction.
    Ny (int): Number of grid points in the y-direction.
    Nz (int): Number of grid points in the z-direction.
    tmax (float): Maximum time for the simulation.
    dt (float): Time step for the simulation.
    N (float): Brunt-Väisälä frequency.
    Umax (float): Maximum velocity.
    latitude (float, optional): Latitude for the Coriolis parameter. Default is 45.
    """
    for setup in setups:
        growth_rate, times, maxV_values, Q, v, U, PSI, zvorticity = eady_analysis(
            Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=latitude, f_cte=setup['f_cte'], 
            perturbation='random', linear_model=setup['linear_model'], rho_cte=setup['rho_cte']
        )
        print(f"growth_rate {setup['label'][-1]}", growth_rate)
        print(f"maxV_values {setup['label'][-1]}", np.max(maxV_values))
        plt.semilogy(times[::20]/3600/24, maxV_values[::20], color=setup['color'], linestyle='--', label=setup['label'])

    plt.semilogy(times[::20]/3600/24, np.exp(0.26*times[::20]/3600/24)/1e2, '-', color='red', label='Analytical Growth Rate')
    plt.grid()
    plt.legend()
    plt.show()
    
setups = [
    {'f_cte': True, 'linear_model': True, 'rho_cte': True, 'color': 'black', 'label': 'Setup 1'},
    {'f_cte': True, 'linear_model': False, 'rho_cte': True, 'color': 'blue', 'label': 'Setup 2'},
    {'f_cte': False, 'linear_model': True, 'rho_cte': True, 'color': 'orange', 'label': 'Setup 3'},
    {'f_cte': False, 'linear_model': False, 'rho_cte': True, 'color': 'grey', 'label': 'Setup 4'},
    {'f_cte': True, 'linear_model': True, 'rho_cte': False, 'color': 'pink', 'label': 'Setup 5'},
    {'f_cte': True, 'linear_model': False, 'rho_cte': False, 'color': 'green', 'label': 'Setup 6'},
    {'f_cte': False, 'linear_model': True, 'rho_cte': False, 'color': 'purple', 'label': 'Setup 7'},
    {'f_cte': False, 'linear_model': False, 'rho_cte': False, 'color': 'yellow', 'label': 'Setup 8'}
]

# # Call the function
# plot_growth_rates(setups, Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45)

########################################################################################

def analyze_and_plot_growth_rates(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45, pert=None):
    """
    Analyze and plot growth rates for different perturbations or random perturbations.

    Parameters:
    Lx, Ly, Lz (float): Domain dimensions.
    Nx, Ny, Nz (int): Number of grid points.
    tmax (float): Maximum simulation time.
    dt (float): Time step.
    N (float): Brunt-Väisälä frequency.
    Umax (float): Maximum velocity.
    latitude (float): Latitude for the Coriolis parameter.
    pert (list or None): List of perturbations or None for random perturbations.
    """
    maxxV_values = []
    if pert is None:
        pert = ['random'] * 20
    else:
        pert = pert + ['random'] * (20 - len(pert))

    for i in range(1, 21):
        U0 = Umax * (z / Lz) ** i
        dU0dz = Umax * i * (z / Lz) ** (i - 1) / Lz
        d2U0dz2 = Umax * i * (i - 1) * (z / Lz) ** (i - 2) / Lz ** 2
        perturbation = pert[i - 1]

        growth_rate, times, maxV_values, Q, v, U, PSI, zvorticity = eady_analysis(
            Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=latitude, f_cte=True, 
            perturbation=perturbation, linear_model=True, rho_cte=True, custom_profile=True, 
            U0=U0, dU0dz=dU0dz, d2U0dz2=d2U0dz2
        )
        maxxV_values.append(np.max(growth_rate))

    plt.semilogy(range(1, 21), maxxV_values, marker='o', color='black', linestyle='--', linewidth=1, markersize=5)
    plt.grid()
    plt.show()

# # Example usage:
# # (Where pert is given by the most unstable perturbations from the previous plot)
# analyze_and_plot_growth_rates(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45, pert=[3, 4, 6, 8, 11, 14, 14, 14, 14, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16])

########################################################################################

def plot_contour_of_v(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45, i=3, perturbation=4):
    """
    Plots the contour of v for a given setup using the Eady model analysis.

    Parameters:
    Lx, Ly, Lz (float): Domain dimensions.
    Nx, Ny, Nz (int): Number of grid points.
    tmax (float): Maximum simulation time.
    dt (float): Time step.
    N (float): Brunt-Väisälä frequency.
    Umax (float): Maximum velocity.
    latitude (float): Latitude for the Coriolis parameter.
    i (int): Exponent for the velocity profile.
    perturbation (int): Perturbation mode.
    """
    z = np.linspace(0, Lz, Nz)
    U0 = Umax * (z / Lz) ** i
    dU0dz = Umax * i * (z / Lz) ** (i - 1) / Lz
    d2U0dz2 = Umax * i * (i - 1) * (z / Lz) ** (i - 2) / Lz ** 2

    growth_rate, times, maxV_values, Q, v, U, PSI, zvorticity = eady_analysis(
        Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=latitude, f_cte=True, 
        perturbation=perturbation, linear_model=True, rho_cte=True, custom_profile=True, 
        U0=U0, dU0dz=dU0dz, d2U0dz2=d2U0dz2
    )

    plt.figure()
    plt.contourf(v[:, -1, :].T, levels=50, cmap='RdBu_r')
    plt.colorbar()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Contour plot of v')
    plt.show()

# Example usage:
# plot_contour_of_v(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45, i=3, perturbation=4)

########################################################################################

def analyze_growth_rate_vs_perturbation(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45, i_range=range(2, 5)):
    """
    Analyze and plot growth rates vs perturbations for different velocity profiles.

    Parameters:
    Lx, Ly, Lz (float): Domain dimensions.
    Nx, Ny, Nz (int): Number of grid points.
    tmax (float): Maximum simulation time.
    dt (float): Time step.
    N (float): Brunt-Väisälä frequency.
    Umax (float): Maximum velocity.
    latitude (float): Latitude for the Coriolis parameter.
    """
    growth_rate_old = 0
    for i in i_range:
        U0 = Umax * (z / Lz) ** i
        dU0dz = Umax * i * (z / Lz) ** (i - 1) / Lz
        d2U0dz2 = Umax * i * (i - 1) * (z / Lz) ** (i - 2) / Lz ** 2
        print('U0 = Umax*(z/Lz)^', i)
        growth_rates = []
        for j in range(1, 16):
            growth_rate, times, maxV_values, Q, v, U, PSI, zvorticity = eady_analysis(
                Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=latitude, f_cte=True, 
                perturbation=j, linear_model=True, rho_cte=True, custom_profile=True, 
                U0=U0, dU0dz=dU0dz, d2U0dz2=d2U0dz2
            )
            
            # Check if the simulation diverged
            if maxV_values[-1] > 1000:
                print('Simulation diverged')
                break
            
            growth_rates.append(growth_rate)
            if growth_rate > growth_rate_old:
                growth_rate_old = growth_rate
                n = j
                
        plt.plot(range(1, len(growth_rates) + 1), growth_rates, marker='o', color='black', linestyle='--', linewidth=1, markersize=5)
        plt.xlabel('Perturbation')
        plt.ylabel('Growth Rate')
        plt.title(f'Growth Rate vs Perturbation (n = {i}) , max growth rate = {growth_rate_old:.3f}')   
        plt.grid(True)
        plt.savefig(f'growth_rate_vs_perturbation_{i}.pdf')

# Example usage:
# analyze_growth_rate_vs_perturbation(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45)

########################################################################################

# PARAMETERS FOR ERA5 DATA
def era5_data_analysis(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=45):
    """
    Perform analysis using ERA5 data and plot the results.

    Parameters:
    Lx, Ly, Lz (float): Domain dimensions.
    Nx, Ny, Nz (int): Number of grid points.
    tmax (float): Maximum simulation time.
    dt (float): Time step.
    N (float): Brunt-Väisälä frequency.
    Umax (float): Maximum velocity.
    latitude (float): Latitude for the Coriolis parameter.
    """
    z = np.linspace(0, Lz, Nz)
    U0 = Umax * (9.43e-01 * (z / Lz) + 5.73e-02 * np.ones_like(z))
    dU0dz = Umax * (9.43e-01 * np.ones_like(z) / Lz)
    d2U0dz2 = Umax * (np.zeros_like(z))

    growth_rate, times, maxV_values, Q, v, U, PSI, zvorticity = eady_analysis(
        Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude=latitude, f_cte=False, 
        perturbation='random', linear_model=False, rho_cte=False, custom_profile=True, U0=U0, dU0dz=dU0dz, d2U0dz2=d2U0dz2
    )

    v_horizontal = np.sqrt(v**2 + U**2)

    lon = np.linspace(-180, 180, Nx)  # Map x to -180° to 180° longitude
    lat = np.linspace(50, 70, Ny)  # Map y to 50° to 70° latitude
    Lon, Lat = np.meshgrid(lon, lat)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8), subplot_kw={'projection': ccrs.Orthographic(central_longitude=0, central_latitude=60)})

    def create_plot(ax, variable=Q, name='Potential Vorticity q (PVU)', arrows=False):
        land = cfeature.NaturalEarthFeature('physical', 'land', '110m', edgecolor='face', facecolor='white')
        ax.add_feature(land, zorder=0)

        contour = ax.contourf(Lon, Lat, variable[:, :, Nz // 2].T, levels=10, cmap='RdBu_r', transform=ccrs.PlateCarree(), zorder=1)

        ax.coastlines(zorder=4)
        gridlines = ax.gridlines(draw_labels=True, zorder=5)

        gridlines.xlabel_style = {'size': 14}
        gridlines.ylabel_style = {'size': 14}

        cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.07, shrink=0.7)  # Increase pad value to move the color bar to the right
        cbar.set_label(name, fontsize=14)
        cbar.ax.tick_params(labelsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)

        if arrows:
            X, Y = np.meshgrid(lon, lat)
            ax.quiver(X, Y, v[:, :, Nz // 2].T, U[:, :, Nz // 2].T, scale=200, color='black', transform=ccrs.PlateCarree())

    create_plot(ax1, name='Potential Vorticity q', arrows=False)
    create_plot(ax2, variable=v_horizontal, name=r'Horizontal velocity perturbation (ms$^-1$)', arrows=False)

    fig.suptitle(f'Time: {times[-1] / 3600 / 24:.0f} days, Max U: {Umax} m/s', fontsize=16, color='black', zorder=6)
    plt.show()

# Example usage:
# era5_data_analysis(Lx=13000000, Ly=2228000, Lz=1e4, Nx=2**6, Ny=2**4, Nz=50, tmax=3600*24*7, dt=800, N=0.01, Umax=50, latitude=45)

########################################################################################
