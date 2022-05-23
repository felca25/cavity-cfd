import numpy as np
from numba import njit

def load_files(Lx, Ly, t, Re):
    
    u = np.loadtxt(f"cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.3f}/3_u.txt")
    v = np.loadtxt(f"cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.3f}/4_v.txt")
    pressure = np.loadtxt(f"cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.3f}/2_pressure.txt")

    return u, v, pressure

@njit
def calculate_stream_function(psi, u, v, N_X, N_Y, D_X, D_Y, TOL):
    
    alpha = -((2.0/(D_X*D_X)) + (2.0/(D_Y*D_Y)))
    
    erro  = 100
    iter = 0

    while erro > TOL:
        R_max = 0
        
        for i in range(1,N_X):
            for j in range(1,N_Y):

                R1 = (-((v[i,j] - v[i-1,j])/D_X) + ((u[i,j] - u[i,j-1])/D_Y))
                
                R2 = (((psi[i+1,j] - 2.0*psi[i,j] + psi[i-1,j])/(D_X*D_X))\
                        +((psi[i,j+1] - 2.0*psi[i,j] + psi[i,j-1])/(D_Y*D_Y)))
                
                R = + R1 - R2
                R = R/alpha
                
                psi[i,j] = psi[i,j] + R
                
        if np.abs(R) > R_max:
            R_max = np.abs(R)
        
        erro = R_max
        iter += 1

    return psi

@njit
def calculate_velocity_plot(u_plot, v_plot, u, v, Nx, Ny):

    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            u_plot[i,j] = 0.5*(u[i,j]+u[i,j-1]) 
            v_plot[i,j] = 0.5*(v[i,j]+v[i-1,j])
    return u_plot, v_plot


        
        