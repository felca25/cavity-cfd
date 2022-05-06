from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import numpy as np



def load_files(Lx, Ly, t, Re):
    
    stream_lines = np.loadtxt(f"cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.3f}/5_stream_function.txt")
    u_plot = np.loadtxt(f"cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.3f}/6_uplot.txt")
    v_plot = np.loadtxt(f"cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.3f}/7_vplot.txt")
    pressure = np.loadtxt(f"cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re}/t_{t:.3f}/2_pressure.txt")

    return u_plot, v_plot, pressure, stream_lines 


if __name__ == "__main__":
    plt.rcParams["font.family"] = "Times New Roman"
    
    Re = 1000
    Lx, Ly = 1., 1.
    t = 2.5
    
    u_plot, v_plot, pressure, stream_lines = load_files(Lx, Ly, t, Re)

    u_plot = u_plot[:-1, :-1]
    v_plot = v_plot[:-1, :-1]
    pressure = pressure[:-2, :-2]
    stream_lines = stream_lines[:-1, :-1]
    
    Nx, Ny = len(u_plot[0]), len(u_plot)
    
    
    velocity = np.zeros([Ny, Nx])
    for i in range(Nx):
        for j in range(Ny):
            velocity[j,i] = np.sqrt(u_plot[j,i]**2 + v_plot[j,i]**2)
            
            
            
    
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    
    Y, X = np.mgrid[0:1:complex(0,Ny), 0:1:complex(0, Nx)]

    print('Plotting graphs ...')
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    title = str(f"Stream Function $\psi$ inside Lid-Driven Cavity with\nRe={Re}, Resolution {Nx}x{Ny} at t={t:.3f}")
    plt.title(title)
    vcmp = ax1.pcolormesh(stream_lines, cmap='jet', rasterized = True, 
                          vmax = stream_lines.max(), vmin = stream_lines.min())
    fig.colorbar(vcmp)
    
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)
    plt.title(f"Pressure inside Lid-Driven Cavity with\nRe={Re}, Resolution {Nx}x{Ny} at t={t:.3f}")
    pssp = ax2.pcolormesh(pressure, cmap='jet', rasterized = True, 
                          vmax = pressure.max(), vmin = pressure.min())
    fig2.colorbar(pssp)
        
    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(111)
    plt.title(f"Streamplot Lid-Driven Cavity with\nRe={Re}, Resolution {Nx}x{Ny} at t={t:.3f}")
    strm = ax3.streamplot(X, Y, u_plot, v_plot, density= 2.5, color= 'k', cmap='jet',
                          linewidth=0.5, arrowsize= 0.5)
    
    fig4 = plt.figure(4)
    ax4 = fig4.add_subplot(111)
    title = str(f"Velocities Modulus inside Lid-Driven Cavity with\nRe={Re}, Resolution {Nx}x{Ny} at t={t:.3f}")
    plt.title(title)
    vcmp = ax4.pcolormesh(velocity, cmap='jet', rasterized = True, 
                          vmax = velocity.max(), vmin = velocity.min())
    fig4.colorbar(vcmp)

    print('Done plotting graphs!')
    plt.show()
        
