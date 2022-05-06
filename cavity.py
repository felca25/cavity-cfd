import numpy as np
import os
import save_data

from functions import *
os.system('cls' if os.name == 'nt' else 'clear')

def run(Nx, Ny, Lx, Ly, reynolds, dx, dy, dt, t_arr, result_params, TOL):


    for i in range(len(reynolds)):
        
        Re = reynolds[i]
        
        save_data.create_paths(Lx, Ly, Re, t_arr, result_params)
        
        # Running the actual simulation from now on
        
        u = np.zeros([Nx+1, Ny+2], float)
        v = np.zeros([Nx+2, Ny+1], float)
        p = np.zeros([Nx+2, Ny+2], float)
        psi = np.zeros([Nx+1,Ny+1])

        u[:,Ny] = 1.0

        u_star = np.copy(u)
        v_star = np.copy(v)

        t = 0
        
        while t <= t_arr[-1]:
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
                if t >= time - TOL and t <= time + TOL:
            
                    rel_path = f'cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.3f}'
            
                    save_data.save_data(aux_arr, result_params, rel_path)
            
        psi = calculate_stream_function(psi, u, v, Nx, Ny, dx, dy, dt, TOL)
