import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

w = 1
Nx = 10
Ny = 10
Y,X = np.mgrid[0:1:complex(0,Ny), 0:1:complex(0, Nx)]

print(X)