import os
import numpy as np
import pandas as pd
os.system('cls' if os.name == 'nt' else 'clear')


def create_paths(Lx, Ly, Re, t_arr, result_params):
    '''
    Creates directories to store cavity flow data
    
    returns:
        all created txt paths
        
    params:
        Re: Reynold's numbers used
        t_arr: time_steps_saved
        result_params: params from the cavity problem 
        to be save
    '''
    paths = []
    
    for j in range(len(t_arr)):
        t_arr[j] = round(t_arr[j], 3)
        path = f'cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t_arr[j]:.3f}'
        
        try:
            os.makedirs(path)
        # os raises an Error when the path already exists, because we don't want our
        # programm to stop, when handle it like so    
        except FileExistsError:
            pass
        
        # Clean-up or Create .txt files to save matrices
        for k in range(len(result_params)):
            
            txt_path = f'{path}/{k}_{result_params[k]}.txt'
            
            f =  open(txt_path, 'w').close()
            paths.append(txt_path)
                
                
    return paths


def save_data(matrix, params, rel_path):
    '''
    Stores numpy array data structures into the disc as .txt files
    
    returns: 1
    
    parameters:
        matrix: list or tuple of the ndarray to be saved
        params: list or tuple with strings of parameters name's to be saved
        rel_path: relative path to the saving destination
    '''
    
    for i in range(len(params)):
        
        path = f'{rel_path}/{i}_{params[i]}.txt'
        
        file = open(path, 'w')
        np.savetxt(path, np.transpose(matrix[i]))
        file.close()

    print(f"Done saving .txt on ... {rel_path}")
    
    return 1


if __name__ == '__main__':
    Nx, Ny = 25, 25
    Lx, Ly = 1., 1.

    Re = (1., 10., 100., 1000.)
    dx, dy = Lx / Nx, Ly / Ny

    for reynolds in Re: 
        dt = 0.25 * reynolds * (dx**2)

    t_final = 0.5
    t_arr = (dt, 10*dt, 25*dt, 100*dt)

    TOL = 1e-5

    result_params = ('u_star', 'v_star', 'pressure', 'u', 'v', 'stream_function')
    
    paths = create_paths(Re, t_arr, result_params)