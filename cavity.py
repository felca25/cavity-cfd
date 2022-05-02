import numpy as np
import os
from matplotlib import pyplot as plt
import save_data
os.system('cls' if os.name == 'nt' else 'clear')

result_params = ('u_star', 'v_star', 'pressure', 'u', 'v', 'stream_function', 'uplot', 'vplot')

def calculate_u_star(u_star, u, v, Nx, Ny, dx, dy, dt, Re):
    
    global se_fudeu
    
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
        
        
    for i in range(1,Nx):
        for j in range(0,Ny):      
            if np.abs(u_star[i,j]) > 1.0:
                print(f'[{i}][{j}]\nu[i,j]={u[i,j]}\nu_star={u_star[i,j]}')
                raise ValueError
            
    return u_star

def calculate_v_star(v_star, v, u, Nx, Ny, dx, dy, dt, Re):
    
    global se_fudeu
    
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
    
    for i in range(0,Nx):
        for j in range(1,Ny):      
            if np.abs(u_star[i,j]) > 1.0:
                print(f'[{i}][{j}]\nu[i,j]={u[i,j]}\nu_star={u_star[i,j]}')
                raise ValueError
        
    return v_star      

def calculate_pressure(p, u_star, v_star, Nx, Ny, dx, dy, dt, TOL):
    global se_fudeu
    error = 100;
    iter = 0
    while error > TOL and iter < 1000:
        # print(iter)
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
                        
                    else:
                        raise IndexError
                    
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
                        
                    else:
                        raise IndexError
                    
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
                    else:
                        raise IndexError
                
                else:
                    raise IndexError
                R = R1 - R2
                R = R/alpha
                p[i,j] = p[i,j] + R
                
                if np.abs(R) > R_max:
                    R_max = np.abs(R)
                    
                if np.abs(p[i,j]) > 1.0:
                    print(f'[{i}][{j}]\np[i,j]={p[i,j]}\n')
                    raise ValueError
                
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

def calculate_new_u(u, u_star, p , Nx, Ny, dx, dt):
    
    for i in range(1,Nx):
        for j in range(-1,Ny+1):
            u[i,j] = u_star[i,j] - dt * ((p[i,j] - p[i-1,j])/dx)
        
    
    for i in range(0,Nx):
        for j in range(0,Ny-1):      
            if np.abs(u_star[i,j]) > 1.0:
                print(f'[{i}][{j}]\nu[i,j]={u[i,j]}\nu_star={u_star[i,j]}')
                raise ValueError

    
    return u

def calculate_new_v(v, v_star, p, Nx, Ny, dx, dt):
    
    for i in range(-1,Nx+1):
        for j in range(1, Ny):
            
            v[i,j] = v_star[i,j] - dt * ((p[i,j] - p[i,j-1])/dy)
            
    for i in range(0,Nx):
        for j in range(0,Ny):      
            if np.abs(u_star[i,j]) > 1.0:
                print(f'[{i}][{j}]\nu[i,j]={u[i,j]}\nu_star={u_star[i,j]}')
                raise ValueError

    
    return v

def calculate_stream_function(psi, u, v, Nx, Ny, dx, dy, dt, TOL):
    
    alpha = -((2.0/(dx*dx)) + (2.0/(dy*dy)))
    
    erro  = 100
    iter = 0
    while erro > TOL:
        R_max = 0
        
        for i in range(1,Nx):
            for j in range(1,Ny):
                R1 = (((v[i,j] - v[i-1,j])/dx) + ((u[i,j] - u[i,j-1])/dy))
                
                R2 = (((psi[i+1,j] - 2.0*psi[i,j] + psi[i-1,j])/(dx*dx))\
                        +((psi[i,j+1] - 2.0*psi[i,j] + psi[i,j-1])/(dy*dy)))
                
                R = - R1 - R2
                R = R/alpha
                
                psi[i,j] = psi[i,j] + R
                
        if np.abs(R) > R_max:
            R_max = np.abs(R)
        
        erro = R_max
        iter += 1    
    return psi

def calculate_velocity_plot(u, v):
    uplot = np.zeros((Nx+1,Ny+1),float) 
    vplot = np.zeros((Nx+1,Ny+1),float) 
    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            uplot[i,j] = 0.5*(u[i,j]+u[i,j-1]) 
            vplot[i,j] = 0.5*(v[i,j]+v[i-1,j])
            
    return uplot, vplot





TOL = 1e-5

Nx, Ny = 25, 25
Lx, Ly = 1., 1.

reynolds = [1, 10, 100, 1000]

dx, dy = Lx / Nx, Ly / Ny

for Re in reynolds:
    dt = 0.25 * Re * (dx**2)

if dt > dx:
    dt *= 0.01

t_arr = (dt, 10*dt, 25*dt, 100*dt, 500*dt)


x = np.linspace(0,Lx,Nx+1)
y = np.linspace(0,Ly,Ny+1)
    
u = np.zeros([Nx+1, Ny+2], float)
v = np.zeros([Nx+2, Ny+1], float)
p = np.zeros([Nx+2, Ny+2], float)
psi = np.zeros([Nx+1,Ny+1])

u[:,Ny] = 1.0

u_star = np.copy(u)
v_star = np.copy(v)

t = 0

save_data.create_paths(Lx, Ly, reynolds, t_arr, result_params)

for k in range(len(reynolds)):
    while t < t_arr[-1]:
        t += dt
        print(t)
        
        u_star = calculate_u_star(u_star, u, v, Nx, Ny, dx, dy, dt, Re)
        
        v_star = calculate_v_star(v_star, v, u, Nx, Ny, dx, dy, dt, Re)
        
        p = calculate_pressure(p, u_star, v_star, Nx, Ny, dx, dy, dt, TOL)
        
        u = calculate_new_u(u, u_star, p, Nx, Ny, dx, dt)
        
        v = calculate_new_v(v, v_star, p, Nx, Ny, dx, dt)
        
        psi = calculate_stream_function(psi, u, v, Nx, Ny, dx, dy, dt, TOL)
        
        uplot, vplot = calculate_velocity_plot(u,v)
        
        aux_arr = (u_star, v_star, p, u, v, psi, uplot, vplot)
        
        
        for time in t_arr:
            if t >= time - 0.1*dt and t <= time + 0.1*dt:
                
                rel_path = f'cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{reynolds[k]}/t_{t:.2f}'
        
                save_data.save_data(aux_arr, result_params, rel_path)
        
psi = calculate_stream_function(psi, u, v, Nx, Ny, dx, dy, dt, TOL)

plt.show()
