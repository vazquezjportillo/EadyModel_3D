import numpy as np 
import matplotlib.pyplot as plt
from functions import transform, inverse_transform, wavenumbers
from scipy.linalg import solve_banded

############################################################################

def eady_model(Lx, Ly, Lz, Nx, Ny, Nz, tmax, dt, N, Umax, latitude = 45, f_cte = True, perturbation = 'random', linear_model = False, rho_cte = True, custom_profile = False, U0 = False, dU0dz = False, d2U0dz2 = False):
    """
    Implements the Eady model for baroclinic instability.
    
    This function simulates the evolution of perturbations in a quasi-geostrophic framework 
    under linear or nonlinear conditions. The Eady model uses spectral methods to compute 
    velocity fields, and vorticity.

    Parameters:
    -----------
    Lx : float
        Domain length in the x-direction (meters).
    Ly : float
        Domain length in the y-direction (meters).
    Lz : float
        Domain height in the z-direction (meters).
    Nx : int
        Number of grid points in the x-direction.
    Ny : int
        Number of grid points in the y-direction.
    Nz : int
        Number of grid points in the z-direction.
    tmax : float
        Maximum simulation time (seconds).
    dt : float
        Time step for the simulation (seconds).
    N : float
        Brunt–Väisälä frequency (rad/s).
    Umax : float
        Maximum zonal wind or flow velocity (m/s).
    latitude : float
        Latitude of the domain (degrees).
    f_cte : bool, optional 
        If True, the Coriolis parameter is constant. Defaults to True. If False, the Coriolis parameter is latitude-dependent and the beta effect is considered.
    perturbation : str or int/float, optional
        Initial perturbation. 'random' for random perturbations, or an integer/float to specify a wave number.
        Defaults to 'random'.
    linear_model : bool, optional
        If True, the model is linear. Defaults to True. If False, the model is nonlinear.
    rho_cte : bool, optional
        If True, the density is constant. Defaults to True. If False, the density is height-dependent.
    custom_profile : bool, optional
        If True, uses a custom vertical velocity profile (U0, dU0dz, d2U0dz2). Defaults to False.
    U0 : array-like, optional
        Custom velocity profile U0(z) if `custom_profile` is True. Defaults to False.
    dU0dz : array-like, optional
        First derivative of the custom velocity profile if `custom_profile` is True. Defaults to False.
    d2U0dz2 : array-like, optional
        Second derivative of the custom velocity profile if `custom_profile` is True. Defaults to False.

    Returns:
    --------
    time : list with time steps
    maxV_values : list with maximum meridional velocity values
    q : ndarray
        Potential vorticity field (real space).
    v : ndarray
        Meridional velocity field (real space).
    u : ndarray
        Zonal velocity field (real space).
    psi : ndarray
        Streamfunction field (real space).
    zeta : ndarray
        Vertical component of vorticity (real space).
    """

    kxx = None
    c = None
    Er = 6371000
    Omega = 7.2921e-5
    g = 9.80665

    def beta(y):
        return 2*Omega/Er*np.cos(np.deg2rad(y))
    
    def f0(y):
        return 2*Omega*np.sin(np.deg2rad(y))
        
    ############################################################################

    # Grid
    x = np.linspace(0,Lx,Nx,endpoint=False)
    y = np.linspace(0,Ly,Ny,endpoint=False)
    z = np.linspace(0,Lz,Nz)
    dz = z[1]-z[0]
    X,Y = np.meshgrid(x,y,indexing='ij')

    # Cte parameter definition
    if f_cte:
        beta = 0
        f0= f0(latitude)
    else:
        beta = beta(latitude)
        f0 = f0(latitude)
    # print('f0:',f0)
        
    # Initial velocity profile
    if custom_profile == False:
        U0 = Umax*(z/Lz)
        dU0dz = Umax/Lz*np.ones_like(z)
        d2U0dz2 = np.zeros_like(z)

    ############################################################################

    kx,ky = wavenumbers(Nx,Ny,Lx,Ly)
    # initial conditions (perturbations only)
    q = np.zeros((Nx,Ny,Nz))
    DPSI_bottom=np.zeros((Nx,Ny))
    DPSI_top=np.zeros((Nx,Ny))
    
    q_hat = transform(q*0)
    if perturbation == 'random':
        
        # Set the random seed for reproducibility
        # seed42 before
        np.random.seed(123)  
                
        DPSI_bottom=np.random.randn(Nx,Ny)
        DPSI_top=np.random.randn(Nx,Ny)
        
    elif type(perturbation) in (int, float): # if perturbation is a number 
        # Forcing just one wave at wavenumber given by perturbation
        kxx=kx[perturbation]
        Ld=N*Lz/f0
        mu=Ld*kxx
        # print('Perturbation wavenumber:',perturbation,'Normalized wavenumber:',mu)
        c=Umax/2+Umax/mu*np.sqrt((mu/2-1/np.tanh(mu/2))*(mu/2-np.tanh(mu/2))+0j) # eady result

        for i in range(Nz):
            q[:,:,i] += 1e-2*(np.cosh(mu*z[i]/Lz)-Umax*c.real/(mu*np.abs(c)**2)*np.sinh(mu*z[i]/Lz))*np.cos(kxx*X) # this is actually the streamfunction
        
        DPSI_bottom=(q[:,:,1]-q[:,:,0])/dz
        DPSI_top=(q[:,:,-1]-q[:,:,-2])/dz
        
    else:
        print('Give a valid perturbation setup. Availables are "random" or an integer')
        exit()

    ############################################################################

    # Fourier space initialization 
    DPSI_bottom_hat = transform(DPSI_bottom)
    DPSI_top_hat=transform(DPSI_top)
    u_hat=np.zeros_like(q_hat)
    v_hat=np.zeros_like(q_hat)
    psi_hatt = np.zeros_like(q_hat)
    zvorticity_hat = np.zeros_like(q_hat)

    ############################################################################

    # Define the banded matrix for the second derivative
    D_banded=np.ones((3,Nz))*(f0/N)**2/dz**2
    D_banded[1]*=-2
    if not rho_cte: # if the density is not constant
        D_banded[0]-=f0**2/(2*dz)/g
        D_banded[2]+=f0**2/(2*dz)/g

    # bottom boundary condition
    D_banded[1,0] = -1/dz
    D_banded[0,1] = 1/dz

    # top boundary condition
    D_banded[1,-1] = 1/dz
    D_banded[2,-2] = -1/dz

    D_banded0=D_banded.copy() # matrix for the zeroth wave number
    D_banded0[2,0]=1 # imposing an extra BC where the streamfunction at the bottom is zero
    D_banded0[1,1]=0
    D_banded0[0,2]=0

    ############################################################################

    def linear_term(q_hat,DPSI_bottom_hat,DPSI_top_hat):
        
        B_hat=np.zeros(Nz,dtype=np.complex128)
        psi_hat=np.zeros_like(B_hat)
        for i in range(Nx//2+1):  # We only need to compute the first half of the spectrum
            for j in range(Ny): 
                if i+j>0:     
                    D_banded[1, 1:-1] = -2*(f0/N)**2/dz**2 - (kx[i]**2 + ky[j]**2)
                    B_hat[1:-1] = q_hat[i,j,1:-1]
                    B_hat[0] = DPSI_bottom_hat[i,j]
                    B_hat[-1] = DPSI_top_hat[i,j]
                    psi_hat[:] = solve_banded((1,1), D_banded, B_hat)
                    v_hat[i,j]=1j*kx[i]*psi_hat
                    np.copyto(psi_hatt[i,j],psi_hat)
                    zvorticity_hat[i,j]=-(kx[i]**2 + ky[j]**2)*psi_hat
                else:
                    B_hat[1:-1] = q_hat[i,j,1:-1]
                    B_hat[0] = DPSI_bottom_hat[i,j]
                    B_hat[-1] = DPSI_top_hat[i,j]
                    B_hat[1]=0
                    v_hat[i,j]=0
                    psi_hatt[i,j]=solve_banded((1,1), D_banded0, B_hat)
                    zvorticity_hat[i,j]=0

        #All terms are linear now 

        s=beta-(f0/N)**2*d2U0dz2
        if not rho_cte:
            s+=f0**2*dU0dz
        nlt_q = U0*(1j*kx*q_hat.T).T + v_hat*s

        
        nlt_bottom = U0[0] * (1j*kx*DPSI_bottom_hat.T).T
        nlt_bottom += v_hat[:,:,0] * (-dU0dz[0])
            
        nlt_top = U0[-1] * (1j*kx*DPSI_top_hat.T).T
        nlt_top += v_hat[:,:,-1] * (-dU0dz[-1])
        V=inverse_transform(v_hat)

        return nlt_q, nlt_bottom, nlt_top, np.max(np.abs(V))

    ############################################################################

    def nonlinear_term(q_hat,DPSI_bottom_hat,DPSI_top_hat):
        
        B_hat=np.zeros(Nz,dtype=np.complex128)
        psi_hat=np.zeros_like(B_hat)
        for i in range(Nx//2+1):  # We only need to compute the first half of the spectrum
            for j in range(Ny): 
                if i+j>0:       
                    D_banded[1, 1:-1] = -2*(f0/N)**2/dz**2 - (kx[i]**2 + ky[j]**2)
                    B_hat[1:-1] = q_hat[i,j,1:-1]
                    B_hat[0] = DPSI_bottom_hat[i,j]
                    B_hat[-1] = DPSI_top_hat[i,j]
                    psi_hat[:] = solve_banded((1,1), D_banded, B_hat)
                    u_hat[i,j]=-1j*ky[j]*psi_hat
                    v_hat[i,j]=1j*kx[i]*psi_hat
                    np.copyto(psi_hatt[i,j],psi_hat)
                    zvorticity_hat[i,j]=-(kx[i]**2 + ky[j]**2)*psi_hat
                else:
                    B_hat[1:-1] = q_hat[i,j,1:-1]
                    B_hat[0] = DPSI_bottom_hat[i,j]
                    B_hat[-1] = DPSI_top_hat[i,j]
                    B_hat[1]=0
                    u_hat[i,j]=0
                    v_hat[i,j]=0
                    psi_hatt[i,j]=solve_banded((1,1), D_banded0, B_hat)
                    zvorticity_hat[i,j]=0

        V=inverse_transform(v_hat,True)

        s=beta-(f0/N)**2*d2U0dz2
        if not rho_cte:
            s+=f0**2*dU0dz
        nlt_q = (U0+inverse_transform(u_hat,True))*inverse_transform((1j*kx*q_hat.T).T,True) 
        nlt_q += V*(inverse_transform(np.moveaxis(1j*ky*np.moveaxis(q_hat,1,-1),-1,1),True) + s )

        nlt_bottom = (U0[0]+inverse_transform(u_hat[:,:,0],True)) * inverse_transform((1j*kx*DPSI_bottom_hat.T).T,True)
        nlt_bottom += V[:,:,0] * (inverse_transform(1j*ky*DPSI_bottom_hat,True)-dU0dz[0])
            
        nlt_top = (U0[-1]+inverse_transform(u_hat[:,:,-1],True)) * inverse_transform((1j*kx*DPSI_top_hat.T).T,True)
        nlt_top += V[:,:,-1] * (inverse_transform(1j*ky*DPSI_top_hat,True)-dU0dz[-1])
                
        return transform(nlt_q,True), transform(nlt_bottom,True), transform(nlt_top,True), np.max(np.abs(V))

    ############################################################################

    def Euler(q_hat,DPSI_bottom_hat,DPSI_top_hat):
        if not linear_model:
            NLT_q_hat,NLT_bottom_hat,NLT_top_hat,maxV=nonlinear_term(q_hat,DPSI_bottom_hat,DPSI_top_hat)
        else:
            NLT_q_hat,NLT_bottom_hat,NLT_top_hat,maxV=linear_term(q_hat,DPSI_bottom_hat,DPSI_top_hat)
        
        q_hat -= dt*NLT_q_hat

        DPSI_bottom_hat -= dt*NLT_bottom_hat
        DPSI_top_hat -= dt*NLT_top_hat
        return maxV

    ############################################################################

    def Leapfrog(q_hat,DPSI_bottom_hat,DPSI_top_hat,q_hat_old,DPSI_bottom_hat_old,DPSI_top_hat_old):
        if not linear_model:
            NLT_q_hat,NLT_bottom_hat,NLT_top_hat,maxV=nonlinear_term(q_hat,DPSI_bottom_hat,DPSI_top_hat)
        else:
            NLT_q_hat,NLT_bottom_hat,NLT_top_hat,maxV=linear_term(q_hat,DPSI_bottom_hat,DPSI_top_hat)
        
        NLT_q_hat[:]= q_hat_old - 2*dt*NLT_q_hat
        
        NLT_bottom_hat[:]= DPSI_bottom_hat_old - 2*dt*NLT_bottom_hat
        NLT_top_hat[:]= DPSI_top_hat_old - 2*dt*NLT_top_hat

        np.copyto(q_hat_old,q_hat)
        np.copyto(DPSI_bottom_hat_old,DPSI_bottom_hat)
        np.copyto(DPSI_top_hat_old,DPSI_top_hat)

        np.copyto(q_hat,NLT_q_hat)
        np.copyto(DPSI_bottom_hat,NLT_bottom_hat)
        np.copyto(DPSI_top_hat,NLT_top_hat)

        return maxV
        
    ############################################################################

    q_hat_old=np.copy(q_hat)
    DPSI_bottom_hat_old=np.copy(DPSI_bottom_hat)
    DPSI_top_hat_old=np.copy(DPSI_top_hat)
    maxV_values = []
    time=[]
    timex=0
    maxV_values.append(Euler(q_hat, DPSI_bottom_hat, DPSI_top_hat))
    timex+=dt;    time.append(timex)
    while timex<tmax:
        maxV_values.append(Leapfrog(q_hat, DPSI_bottom_hat, DPSI_top_hat, q_hat_old, DPSI_bottom_hat_old, DPSI_top_hat_old))
        timex+=dt;    time.append(timex)
        if np.max(maxV_values)>10:
            print('Simulation diverged');            break
            

    # YY=np.repeat(y[:,np.newaxis],Nz,axis=1)/(N*Lz/f0)
    YY=np.repeat(y[:,np.newaxis],Nz,axis=1)
    # YY = np.where(YY == 0, 1, YY)
    # print(YY)
    Q0=f0+(beta-((f0/N)**2)*d2U0dz2)*YY
    PSI0=-YY*U0
    if not rho_cte:
        Q0+=f0**2*dU0dz*YY/g
    # YY=np.repeat(y[:,np.newaxis],Nz,axis=1)
    # Q0=(f0+((f0/N)**2)*d2U0dz2)/YY
    # # print(f0*10**6)
    # PSI0=-U0/YY

    # Q0 = 0
    
    return time, maxV_values, inverse_transform(q_hat)+Q0, inverse_transform(v_hat), inverse_transform(u_hat)+U0, inverse_transform(psi_hatt)+PSI0, inverse_transform(zvorticity_hat)