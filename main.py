import numpy as np
import os
from functions import *
import cavity
os.system('cls' if os.name == 'nt' else 'clear')

result_params = ('u_star', 'v_star', 'pressure', 'u', 'v', 'stream_function', 'uplot', 'vplot')

TOL = 1e-5

Nx, Ny = 25, 25
Lx, Ly = 1., 1.

reynolds = [1, 10, 100, 1000]

dx, dy = Lx / Nx, Ly / Ny

t_arr_mult = (1., 10., 25., 100., 500.)

cavity.run(Nx, Ny, Lx, Ly, reynolds, dx, dy, t_arr_mult, result_params)