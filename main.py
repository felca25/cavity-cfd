import numpy as np

import os
from functions import *
import cavity
os.system('cls' if os.name == 'nt' else 'clear')

result_params = ('u_star', 'v_star', 'pressure', 'u', 'v', 'stream_function', 'uplot', 'vplot')

TOL = 1e-8

Nx, Ny = 100, 100
Lx, Ly = 1., 1.

reynolds = (10,)

dx, dy = Lx / Nx, Ly / Ny

dt = 0.001

# we need to first apply the stability and
# convergence conditions so that the answers don't blow up
# this should be at the main function because we ought to want
# to have fine manual control over the process. But having this here
# Ensures that the values won't blow up as those variables are Re dependent
# What will be defined are the multipliers of the time steps to be saved

if dx > 1 / np.sqrt(reynolds[0]):
    dx = (1 / np.sqrt(reynolds[0])) - TOL
    
if dy > 1 / np.sqrt(reynolds[0]):
    dy = (1 / np.sqrt(reynolds[0])) - TOL

if dt > 0.25 * reynolds[0] * (dx**2) or dt > 0.25 * reynolds[0] * (dy**2):
    dl = min([dx, dy])
    dt = 0.25 * reynolds[0] * (dl**2)

if dt > dx:
    dt = dx - TOL


t_mult = [1., 5., 10., 25., 50., 75., 100., 250., 500., 750., 1000.,
          1250., 1500., 1750., 2000., 2500., 5000., 7500., 10000., 12500,
          15000, 17500]

t_arr = np.zeros(len(t_mult))
for i in range(len(t_mult)):
    t_arr[i] = t_mult[i] * dt
    

cavity.run(Nx, Ny, Lx, Ly, reynolds, dx, dy, dt, t_arr, result_params, TOL)