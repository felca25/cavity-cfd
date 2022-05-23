import numpy as np
import os
import save_data

from functions import *
os.system('cls' if os.name == 'nt' else 'clear')

def run(N_X, N_Y, L_X, L_Y, REYNOLDS, D_X, D_Y, D_T, T_MULT, RESULT_PARAMETERS, TOL):
    # we need to first apply the stability and
    # convergence conditions so that the answers don't blow up
    # this should be at the main function because we ought to want
    # to have fine manual control over the process. But having this here
    # Ensures that the values won't blow up as those variables are Re dependent
    # What will be defined are the multipliers of the time steps to be saved
    convergence_check(REYNOLDS, D_X, D_Y, D_T, TOL)
    
    t_arr = create_t_arr(D_T, T_MULT)
    
    save_data.create_paths(t_arr, L_X, L_Y, REYNOLDS, RESULT_PARAMETERS)
    
    # Running the actual simulation from now on
    
    u = np.zeros([N_X+1, N_Y+2], float)
    v = np.zeros([N_X+2, N_Y+1], float)
    p = np.zeros([N_X+2, N_Y+2], float)

    u[:,N_Y] = 1.0

    u_star = np.copy(u)
    v_star = np.copy(v)

    t = 0
    
    while t <= t_arr[-1]:
        t += D_T
        print(t)
        
        u_star = calculate_u_star(u_star, u, v, N_X, N_Y, D_X, D_Y, D_T, REYNOLDS)
        
        v_star = calculate_v_star(v_star, v, u, N_X, N_Y, D_X, D_Y, D_T, REYNOLDS)
        
        p = calculate_pressure(p, u_star, v_star, N_X, N_Y, D_X, D_Y, D_T, TOL)
        
        u = calculate_new_u(u, u_star, p, N_X, N_Y, D_X, D_T)
        
        v = calculate_new_v(v, v_star, p, N_X, N_Y, D_X, D_T)
        
        aux_arr = (u_star, v_star, p, u, v)
        
        
        for time in t_arr:
            if t >= time - TOL and t <= time + TOL:
        
                rel_path = f'cavity_results/{L_X:.2f}x{L_Y:.2f}/Re_{REYNOLDS}/t_{t:.3f}'
        
                save_data.save_data(aux_arr, RESULT_PARAMETERS, rel_path)
        

