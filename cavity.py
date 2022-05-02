import numpy as np
import os
import save_data

from functions import *
os.system('cls' if os.name == 'nt' else 'clear')

def run(Nx, Ny, Lx, Ly, reynolds, dx, dy, t_arr_mult, result_params):


    for i in range(len(reynolds)):
        
        Re = reynolds[i]
        
        if dx > 1 / np.sqrt(Re):
            dx = (1 / np.sqrt(Re)) - TOL
            
        dt = 0.25 * Re * (dx**2)
        
        if dt > dx:
            dt *= 0.01
            
        t_arr = np.zeros(len(t_arr_mult))
        
        for i in range(len(t_arr)):
            t_arr[i] = t_arr_mult[i] * dt
        
        save_data.create_paths(Lx, Ly, Re, t_arr, result_params)
        
        u = np.zeros([Nx+1, Ny+2], float)
        v = np.zeros([Nx+2, Ny+1], float)
        p = np.zeros([Nx+2, Ny+2], float)
        psi = np.zeros([Nx+1,Ny+1])

        u[:,Ny] = 1.0

        u_star = np.copy(u)
        v_star = np.copy(v)

        t = 0
        
        while t < t_arr[-1]:
            t += dt
            print(t)
            
            u_star = calculate_u_star(u_star, u, v, Nx, Ny, dx, dy, dt, Re)
            
            v_star = calculate_v_star(v_star, v, u, Nx, Ny, dx, dy, dt, Re)
            
            p = calculate_pressure(p, u_star, v_star, Nx, Ny, dx, dy, dt, TOL)
            
            u = calculate_new_u(u, u_star, p, Nx, Ny, dx, dt)
            
            v = calculate_new_v(v, v_star, p, Nx, Ny, dx, dt)
            
            psi = calculate_stream_function(psi, u, v, Nx, Ny, dx, dy, dt, TOL)
            
            uplot, vplot = calculate_velocity_plot(u, v, Nx, Ny)
            
            aux_arr = (u_star, v_star, p, u, v, psi, uplot, vplot)
            
            
            for time in t_arr:
                if t >= time - 0.1*dt and t <= time + 0.1*dt:
                    
                    rel_path = f'cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.2f}'
            
                    save_data.save_data(aux_arr, result_params, rel_path)
            
        psi = calculate_stream_function(psi, u, v, Nx, Ny, dx, dy, dt, TOL)

if __name__ == "__main__":
    
    result_params = ('u_star', 'v_star', 'pressure', 'u', 'v', 'stream_function', 'uplot', 'vplot')

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
    
    run(Nx, Ny, Lx, Ly, reynolds, dx, dy, dt, t_arr, result_params)