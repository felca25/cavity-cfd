import numpy as np
import os
import time as tm
import save_data

from functions import *
os.system('cls' if os.name == 'nt' else 'clear')

def run(INITIAL_PARAMETERS, save_results=True):
    """_summary_
    Solves and saves the results for the lid-driven cavity problem for the initial arguments given.

    Args:
        INITIAL_PARAMETERS (array like): Here the parameters should be passed in the following order:
        
            N_X, N_Y, L_X, L_Y, REYNOLDS, D_X, D_Y, D_T, t_mult, RESULT_PARAMETERS, TOL
        
        for:
            N_X (int): number of cells in the x axis
            N_Y (int): number of cells in the y axis
            L_X (float): characteristic length of the x axis
            L_Y (float): characteristic length of the y axis
            REYNOLDS (tuple): Reynold's numbers to solve the system for
            D_X (float): spacial x step between discretized points
            D_Y (float): spacial y step between discretized points
            D_Y (float): time step for a given Reynold's number
            t_mult (tuple): tuple of factors to multiply dt to save the results
            RESULT_PARAMETERS (tuple): tuple of which parameters to save under the specified name
            TOL (float): tolerance of the solver (should be a value smaller 1e-5 and larger than 1e-9)
    
    KArgs:
        save_results (bool): saves results by default, set to False if saving is not desirable

    Raises:
        IndexError: when Reynold's numbers tuple and time steps tuples don't match in shape

    Returns:
        boolean : returns True whenever the functions has ran propperly
    """
    
    N_X, N_Y, L_X, L_Y, REYNOLDS,\
        D_X, D_Y, D_T, T_MULT, RESULT_PARAMETERS, TOL = INITIAL_PARAMETERS
        
    if len(REYNOLDS) > len(D_T):
        raise IndexError('Time step array and Reynolds array must be the same size')
    
    save_data.save_initial_conditions(INITIAL_PARAMETERS)
    
    time_taken = []
    for i, Re in enumerate(REYNOLDS):
        dt = D_T[i]
        start_time = tm.time()
        
        # we need to first apply the stability and
        # convergence conditions so that the answers don't blow up
        # this should be at the main function because we ought to want
        # to have fine manual control over the process. But having this here
        # Ensures that the values won't blow up as those variables are Re dependent
        # What will be defined are the multipliers of the time steps to be saved
            
        convergence_check(Re, D_X, D_Y, dt, TOL)
        
        t_arr = create_t_arr(dt, T_MULT)
        
        save_data.create_paths(t_arr, L_X, L_Y, Re, RESULT_PARAMETERS)
        
        # Running the actual simulation from now on
        
        u = np.zeros([N_X+1, N_Y+2], float)
        v = np.zeros([N_X+2, N_Y+1], float)
        p = np.zeros([N_X+2, N_Y+2], float)

        u[:,N_Y] = 1.0

        u_star = np.copy(u)
        v_star = np.copy(v)

        t = 0
        
        while t <= t_arr[-1]:
            t += dt
            print(t)
            
            u_star = calculate_u_star(u_star, u, v, N_X, N_Y, D_X, D_Y, dt, Re)
            
            v_star = calculate_v_star(v_star, v, u, N_X, N_Y, D_X, D_Y, dt, Re)
            
            p = calculate_pressure(p, u_star, v_star, N_X, N_Y, D_X, D_Y, dt, TOL)
            
            u = calculate_new_u(u, u_star, p, N_X, N_Y, D_X, dt)
            
            v = calculate_new_v(v, v_star, p, N_X, N_Y, D_X, dt)
            
            aux_arr = (u_star, v_star, p, u, v)
            
            if save_results:
                for time in t_arr:
                    if t >= time - TOL and t <= time + TOL:
                
                        rel_path = f'cavity_results/{L_X:.2f}x{L_Y:.2f}/Re_{Re}/t_{t:.3f}'
                
                        save_data.save_data(aux_arr, RESULT_PARAMETERS, rel_path)
        
        end_time = tm.time()
        time_taken.append(end_time - start_time)
        
        print(f'Done solving for Reynolds {Re}')
        print(f'Took {time_taken[i]} s')
    
    print(f'Done solving system')
    
    return True
        

