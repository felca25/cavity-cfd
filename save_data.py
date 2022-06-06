import os
import json
import numpy as np

os.system('cls' if os.name == 'nt' else 'clear')


def create_paths(t_arr: list, Lx: float, Ly: float, Re: float, result_params: tuple, clean_up=True):
    """_summary_
        Creates directorys and files to store the results from the cavity problem
    Args:
        t_arr (list): time_steps_saved
        Lx (float): _description_
        Ly (float): _description_
        Re (float): Reynold's numbers used
        result_params (tuple): parameters to be saved and their respective names
        clean_up (bool, optional): controls the option to clean-up the previous stored data from files. Defaults to True.

    Returns:
        list: with all pathnames created for the .txt files 
    """

    path_list = []

    for time_saved in t_arr:
        time_saved = round(time_saved, 3)
        path = f'cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{time_saved:.3f}'

        try:
            os.makedirs(path)
        # os raises an Error when the path already exists, because we don't want our
        # programm to stop, when handle it like so
        except FileExistsError:
            pass

        # Clean-up or Create .txt files to save matrices
        if clean_up:
            for k, result_param in enumerate(result_params):

                txt_path = f'{path}/{k}_{result_param}.txt'

                open(txt_path, 'w', encoding='utf-8').close()
                path_list.append(txt_path)

    return path_list


def save_data(matrix: np.ndarray, params: tuple, rel_path: str):
    """_summary_
    Saves all the parameters given in the matrix into their respective .txt
    files using np.savetxt

    Args:
        matrix (tuple): tuple with ndarray's results to be saved
        params (tuple): tuple with parameters to be saved
        rel_path (str): relative path

    Returns:
        bool: returns True if saving process has gone right
    """

    for i, parameter in enumerate(params):

        path = f'{rel_path}/{i}_{parameter}.txt'

        file = open(path, 'w', encoding='utf-8')
        np.savetxt(path, np.transpose(matrix[i]))
        file.close()

    print(f"Done saving .txt on ... {rel_path}")

    return True


def save_initial_conditions(initial_parameters):
    """_summary_

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

    Returns:
        bool: returns True if the process has gone right
    """

    N_X, N_Y, L_X, L_Y, REYNOLDS, D_X, D_Y, time_steps,\
        t_mult, RESULT_PARAMETERS, TOL = initial_parameters

    for i, Re in enumerate(REYNOLDS):
        parameters = {
            'N_X': N_X,
            'N_Y': N_Y,
            'L_X': L_X,
            'L_Y': L_Y,
            "Reynold's number": Re,
            'D_X': D_X,
            'D_Y': D_Y,
            'time_step': time_steps[i],
            't_mult': t_mult,
            'RESULT_PARAMETERS': RESULT_PARAMETERS,
            'TOLERANCE': TOL,
        }
        data = {'initial conditions': parameters}

        with open(f'cavity_results/{L_X:.2f}x{L_Y:.2f}/Re_{int(Re)}/initial_conditions.json', 'w', encoding='ascii') as outfile:
            json.dump(data, outfile, indent=4)
            print(
                f'Successfuly saved initial conditions for Re {Re} data in initial_conditions.json file')

    return True


if __name__ == '__main__':
    from main import INITIAL_PARAMETERS
    save_initial_conditions(INITIAL_PARAMETERS)
