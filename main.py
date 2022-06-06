import os
import functions.cavity as cavity
os.system('cls' if os.name == 'nt' else 'clear')

RESULT_PARAMETERS = ('u_star', 'v_star', 'pressure', 'u', 'v')

TOL = 1e-8

N_X, N_Y = 100, 100
L_X, L_Y = 1., 1.

REYNOLDS = (1000,)

D_X, D_Y = L_X / N_X, L_Y / N_Y

D_T = (0.01,)

t_mult = [1., 5., 10., 25., 50., 75., 100., 250., 500., 750., 1000.,
        1250., 1500., 1750., 2000., 2500., 5000., 7500., 10000., 12500,
        15000, 17500]

INITIAL_PARAMETERS = (N_X, N_Y, L_X, L_Y, REYNOLDS, D_X, D_Y, D_T, t_mult, RESULT_PARAMETERS, TOL)

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

if __name__ == "__main__":
    cavity.run(INITIAL_PARAMETERS)
