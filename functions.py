import numpy as np
from numba import njit

@njit
def calculate_u_star(u_star, u, v, Nx, Ny, dx, dy, dt, Re):
    """_summary_
    
    Calculates the intermediary velocity field u_star

    Args:
        u_star (class 'numpy.ndarray'): intermediary velocity field u* in the x axis
        u (class 'numpy.ndarray'): velocity field u
        v (class 'numpy.ndarray'): velocity field v
        Nx (int): numbers of cells in the x axis
        Ny (int): number of cells in the y axis
        dx (float): step distance in the x axis
        dy (float): step distance in the y axis
        dt (float): step in time
        Re (float): Reynold's number

    Returns:
        _type_: class 'numpy.ndarray'
    """
    
    for i in range(1,Nx):
        for j in range(0,Ny):        
            C1 = 0.25 * (v[i, j+1] + v[i-1, j+1] + v[i,j] + v[i-1,j])

            R_conv = dt * ((u[i,j] * ((u[i+1,j] - u[i-1,j]) / (2*dx)))
                        + (C1 * ((u[i,j+1] - u[i,j-1]) / (2*dy))))
            
            R_dif = (dt/Re) * (((u[i+1,j] - 2*u[i,j] + u[i-1,j]) / (dx*dx))
                                + ((u[i,j+1] - 2*u[i,j] + u[i,j-1]) / (dy*dy)))
            
            u_star[i,j] = u[i,j] - R_conv + R_dif

    for g in range(0, Ny):
        u_star[0,g] = 0.0
        u_star[Nx,g] = 0.0
        
    for k in range(0, Nx+1):
        u_star[k,-1] = - u_star[k,0]
        u_star[k,Ny] = 2.0 - u_star[k, Ny-1]
    
             
    return u_star

@njit
def calculate_v_star(v_star, v, u, Nx, Ny, dx, dy, dt, Re):
    
    for i in range(0,Nx):
        for j in range(1,Ny):     
            
            C2 = 0.25 * (u[i+1, j] + u[i, j] + u[i+1,j-1] + u[i,j-1])

            R_conv = -dt * (C2 * ((v[i+1,j] - v[i-1,j]) / (2*dx)))\
                        -dt * (v[i,j] * ((v[i,j+1] - v[i,j-1]) / (2*dy)))
            
            R_dif = (dt/Re) * ((v[i+1,j] - 2*v[i,j] + v[i-1,j]) / (dx*dx))\
                       +(dt/Re) * ((v[i,j+1] - 2*v[i,j] + v[i,j-1]) / (dy*dy))
            
            v_star[i,j] = v[i,j] + R_conv + R_dif
            
    for j in range(0, Ny+1):
        v_star[-1,j] = - v_star[0,j]
        v_star[Nx,j] = - v_star[Nx-1,j]
    
    for i in range(0, Nx):
        v_star[i,0] = 0.0
        v_star[i,Ny] = 0.0

        
    return v_star     
 
@njit
def calculate_pressure(p, u_star, v_star, Nx, Ny, dx, dy, dt, TOL):
    error = 100
    iter = 0
    while error > TOL and iter < 1000:
        R_max = 0
        for i in range(0,Nx):
            for j in range(0,Ny):
                       
                R1 = ((u_star[i+1, j] - u_star[i,j]) / (dt*dx))\
                        +((v_star[i, j+1] - v_star[i,j]) / (dt*dy))
                    
                if i == 0:
                    if  j == 0:
                        
                        alpha  = - ((1.0/(dx*dx)) + (1.0/(dy*dy)))
                        
                        R2 = (((p[i+1,j] - 1.0*p[i,j]) / (dx*dx))
                                +((p[i,j+1] - 1.0*p[i,j]) / (dy*dy)))

                    elif j == Ny-1:
                        
                        alpha  = - ((1.0/(dx*dx)) + (1.0/(dy*dy)))         
                        
                        R2 = (((p[i+1,j] - 1.0*p[i,j]) / (dx*dx))
                                +((-1.0*p[i,j] + p[i,j-1]) / (dy*dy)))
                        
                    elif j != 0  and j != Ny-1:
                        
                        alpha = -((1.0/(dx*dx)) + (2.0/(dy*dy)))         
                    
                        R2 = (((p[i+1,j] - 1.0*p[i,j]) / (dx*dx))
                                +((p[i,j+1] - 2.0*p[i,j] + p[i,j-1]) / (dy*dy)))
                    
                elif i == Nx-1:
                    if  j == 0:
                        
                        alpha  = - ((1.0/(dx*dx)) + (1.0/(dy*dy)))
                        
                        R2 = (((-1.0*p[i,j] + p[i-1,j]) / (dx*dx))
                                +((p[i,j+1] - 1.0*p[i,j]) / (dy*dy)))
                        
                    elif j == Ny-1:
                        
                        alpha  = - ((1.0/(dx*dx)) + (1.0/(dy*dy)))         
                        
                        R2 = (((-1.0*p[i,j] + p[i-1,j]) / (dx*dx))
                                +((-1.0*p[i,j] + p[i,j-1]) / (dy*dy)))
                        
                    elif j != 0  and j != Ny-1:
                        
                        alpha = -((1.0/(dx*dx)) + (2.0/(dy*dy)))         
                    
                        R2 = (((-1.0*p[i,j] + p[i-1,j]) / (dx*dx))
                                +((p[i,j+1] - 2.0*p[i,j] + p[i,j-1]) / (dy*dy)))
                    
                elif i != 0 and i != Nx-1:
                    if  j == 0:
            
                        alpha  = - ((2.0/(dx*dx)) + (1.0/(dy*dy)))
                        
                        R2 = (((p[i+1,j] - 2.0*p[i,j] + p[i-1,j]) / (dx*dx))
                                +((p[i,j+1] - 1.0*p[i,j]) / (dy*dy)))
                
                    elif j == Ny-1:
                        
                        alpha  = - ((2.0/(dx*dx)) + (1.0/(dy*dy)))         
                        
                        R2 = (((p[i+1,j] - 2.0*p[i,j] + p[i-1,j]) / (dx*dx))
                                +((-1.0*p[i,j] + p[i,j-1]) / (dy*dy)))
                    
                    elif j != 0  and j != Ny-1:
                        
                        alpha = -((2.0/(dx*dx)) + (2.0/(dy*dy)))         
                    
                        R2 = (((p[i+1,j] - 2.0*p[i,j] + p[i-1,j]) / (dx*dx))
                                +((p[i,j+1] - 2.0*p[i,j] + p[i,j-1]) / (dy*dy)))
                
                
                R = R1 - R2
                R = R/alpha
                p[i,j] = p[i,j] + R
                
                if np.abs(R) > R_max:
                    R_max = np.abs(R)
                
        for i in range(0,Nx):
            p[i,-1] = p[i,0]
            p[i,Ny] = p[i,Ny-1]
            
        for j in range(0,Ny):
            p[-1,j] = p[0,j]
            p[Nx,j] = p[Nx-1,j]
        
        p[-1,-1] = p[0,0]
        p[-1,Ny] = p[0,Ny-1]
        p[Nx,-1] = p[Nx-1,0]
        p[Nx,Ny] = p[Nx-1,j]
                    
        error = R_max
        iter += 1
            
    return p

@njit
def calculate_new_u(u, u_star, p , Nx, Ny, dx, dt):
    
    for i in range(1,Nx):
        for j in range(-1,Ny+1):
            
            u[i,j] = u_star[i,j] - dt * ((p[i,j] - p[i-1,j])/dx)

    return u
@njit
def calculate_new_v(v, v_star, p, Nx, Ny, dy, dt):
    
    for i in range(-1,Nx+1):
        for j in range(1, Ny):
            
            v[i,j] = v_star[i,j] - dt * ((p[i,j] - p[i,j-1])/dy)
            
    return v

def convergence_check(Re, dx, dy, dt, TOL):
    
    reynolds_inv_sqrt =  1 / np.sqrt(Re)
    dl = min([dx, dy])   

    if dx > reynolds_inv_sqrt:
        raise ValueError(f"For convergence the value of dx must be smaller than {reynolds_inv_sqrt}")
        
    if dy > reynolds_inv_sqrt:
        raise ValueError(f"For convergence the value of dy must be smaller than {reynolds_inv_sqrt}")

    if dt > 0.25 * Re * (dl**2):
        dt = 0.25 * Re * (dl**2)
        raise ValueError(f"For convergence the value of dt must be smaller than {dt}")

    if dt > dx:
        dt = dx - TOL
        raise ValueError(f"For convergence the value of dt must be smaller than {dx}")
    
    return 1



def create_t_arr(dt, t_mult):
    t_arr = np.zeros(len(t_mult))
    
    for i, factor in enumerate(t_mult):
        t_arr[i] = factor * dt

    return t_arr
