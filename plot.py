from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import numpy as np
from functions.plotting_functions import *

plt.rcParams["font.family"] = "Times New Roman"

def plot_cavity(time, u_plot, v_plot, pressure, psi, velocity, Re, N_X, N_Y, L_X, L_Y):
            
    
    x = np.linspace(0, L_X, N_X)
    y = np.linspace(0, L_Y, N_Y)
    
    Y, X = np.mgrid[0:1:complex(0, N_Y), 0:1:complex(0, N_X)]

    print('Plotting graphs ...')
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    title = str(f"Re={Re} at t={time:.3f}")
    ax1.set_title(title)
    vcmp = ax1.contourf(X, Y, psi, levels= 10, cmap='jet', 
                        vmax = psi.max(), vmin = psi.min(), linestyles="dashed")
    vcmp2 = ax1.contour(X, Y, psi, levels=10, colors='k')
    ax1.clabel(vcmp2, fmt='%2.2f', colors='k', fontsize=10)
    ax1.set_xlabel('$x$ axis')
    ax1.set_ylabel('$y$ axis')
    fig.colorbar(vcmp)
    plt.savefig(f"Results_{L_X:.2f}x{L_Y:.2f}/StreamFunction_Re={Re}_{N_X}x{N_Y}_t={time:.3f}.png")
    
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)
    ax2.set_title(f"Re={Re} at t={time:.3f}")
    pssp = ax2.pcolormesh(pressure, cmap='jet', rasterized = True, 
                        vmax = pressure.max(), vmin = pressure.min())
    ax2.set_xlabel('Grid point at $x$ axis')
    ax2.set_ylabel('Grid point at $y$ axis')
    fig2.colorbar(pssp)
    plt.savefig(f"Results_{L_X:.2f}x{L_Y:.2f}/Pressure_Re={Re}_{N_X}x{N_Y}_t={time:.3f}.png")
        
    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(111)
    ax3.set_title(f"Re={Re} at t={time:.3f}")
    strm = ax3.streamplot(X, Y, u_plot, v_plot, density= 2.5, color= 'k', cmap='jet',
                        linewidth=0.5, arrowsize= 0.5)
    ax3.set_xlabel('$x$ axis')
    ax3.set_ylabel('$y$ axis')
    plt.savefig(f"Results_{L_X:.2f}x{L_Y:.2f}/StreamPlot_Re={Re}_{N_X}x{N_Y}_t={time:.3f}.png")
    
    fig4 = plt.figure(4)
    ax4 = fig4.add_subplot(111)
    title = str(f"Re={Re} at t={time:.3f}")
    ax4.set_title(title)
    vcmp = ax4.pcolormesh(velocity, cmap='jet', rasterized = True, 
                        vmax = velocity.max(), vmin = velocity.min())
    ax4.set_xlabel('Grid point at $x$ axis')
    ax4.set_ylabel('Grid point at $y$ axis')
    fig4.colorbar(vcmp)
    plt.savefig(f"Results_{L_X:.2f}x{L_Y:.2f}/Velocity_Modulus_Re={Re}_{N_X}x{N_Y}_t={time:.3f}.png")

    print('Done plotting graphs!')
    
    plt.close('all')
    
def calculate_line_plot(u_plot, v_plot, N_X, N_Y, L_X, L_Y, D_X, D_Y):
    x_array = np.arange(0, L_X, D_X)
    y_array = np.arange(0, L_Y, D_Y)
    
    u_center_line = np.zeros(len(y_array))
    v_center_line = np.zeros(len(x_array))

    x_position = N_X//2
    y_position = N_Y//2
    
    for j in range(N_Y):
        u_center_line[j] = u_plot[x_position, j]
    for i in range(N_X):
        v_center_line[i] = v_plot[i, y_position]
    
    return x_array, y_array, u_center_line, v_center_line
    
def plot_velocity_line(x_array, y_array, u_center_line, v_center_line, time, N_X, N_Y, L_X, L_Y, Re):
    x_position = N_X//2
    y_position = N_Y//2
    
    fig = plt.figure('U by Y Line Plot')
    ax = fig.add_subplot(111)
    ax.plot(y_array,u_center_line )
    title = str(f"$u$ at $y$={x_array[x_position]}")
    ax.set_title(title)
    ax.set_xlabel('Velocity $u$')
    ax.set_ylabel('$y$ axis')
    plt.savefig(f"Results_{L_X:.2f}x{L_Y:.2f}/Velocity_U_Line_Re={Re}_{N_X}x{N_Y}_t={time:.3f}.png")
    
    fig2 = plt.figure('X by V Line Plot')
    ax2 = fig2.add_subplot(111)
    ax2.plot(x_array, v_center_line)
    title = str(f"$v$ at $x$={y_array[y_position]}")
    ax2.set_title(title)
    ax2.set_xlabel('$x$ axis')
    ax2.set_ylabel('Velocity $v$')
    plt.savefig(f"Results_{L_X:.2f}x{L_Y:.2f}/Velocity_V_Line_Re={Re}_{N_X}x{N_Y}_t={time:.3f}.png")
    
    plt.close('all')
    
def calculate_velocity_plot(u_plot, v_plot, u, v, Nx, Ny):
    
    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            u_plot[i,j] = 0.5*(u[i,j]+u[i,j-1]) 
            v_plot[i,j] = 0.5*(v[i,j]+v[i-1,j])
    return u_plot, v_plot

def calculate_velocity_modulus(velocity, u_plot, v_plot, N_X, N_Y):
    
    for i in range(N_X):
        for j in range(N_Y):
            velocity[i,j] = np.sqrt(u_plot[i,j]**2 + v_plot[i,j]**2)       
        


if __name__ == "__main__":
    
    from main import INITIAL_PARAMETERS
    
    D_T = INITIAL_PARAMETERS[8]
    TOL = INITIAL_PARAMETERS[-1]
    
    L_X, L_Y = 1., 1.
    N_X, N_Y = 100, 100
    
    D_X, D_Y = L_X / N_X, L_Y / N_Y

    REYNOLDS = (1000,)

    times = [50.,]
    # 0.438,
    
    for j, Re in enumerate(REYNOLDS):
        u, v, pressure = load_files(L_X, L_Y, times[j], Re)
        u = np.transpose(u)
        v = np.transpose(v)
        

        u_plot = np.zeros([N_X+1,N_Y+1],float) 
        v_plot = np.zeros([N_X+1,N_Y+1],float)
        psi = np.zeros([N_X+1,N_Y+1], float)
        velocity = np.zeros([N_X, N_Y])
        
        
        print("Calculating Stream Function ...")
        psi = calculate_stream_function(psi, u, v, N_X, N_Y, D_X, D_Y, TOL)
        print("Done Calculating Stream Function !!!") 
            

        print('Calculating Velocity Values for Plotting ...')
        u_plot, v_plot = calculate_velocity_plot(u_plot, v_plot, u, v, N_X, N_Y)
        print("Done Calculating Velocity Values to plot !!!")
        
        u_plot = u_plot[:-1, :-1]
        u = u[:-1, :-2]
        v_plot = v_plot[:-1, :-1]
        pressure = pressure[:-2, :-2]
        psi = psi[:-1, :-1]
        
        calculate_velocity_modulus(velocity, u_plot, v_plot, N_X, N_Y)
        
        print(N_X, N_Y)
        print(len(u[0]), len(u))
        print(len(v[0]), len(v))
        print(len(psi[0]), len(psi))
        print(len(u_plot[0]), len(u_plot))
        print(len(v_plot[0]), len(v_plot))
        
        x_array, y_array, u_center_line, v_center_line = calculate_line_plot(u_plot, v_plot, N_X, N_Y, L_X, L_Y, D_X, D_Y)
        plot_velocity_line(x_array, y_array, u_center_line, v_center_line, times[j], N_X, N_Y, L_X, L_Y, Re)
        

        u_plot = np.transpose(u_plot)
        v_plot = np.transpose(v_plot)
        psi = np.transpose(psi)
        velocity = np.transpose(velocity)
        
        plot_cavity(times[j], u_plot, v_plot, pressure, psi, velocity, Re, N_X, N_Y, L_X, L_Y)

        
