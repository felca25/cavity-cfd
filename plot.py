from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np


def read_np_file(filepath):
    file = open(filepath, 'r')

    lines = []
    for line in file:
        lines.append(line)
        
    for i in range(0, len(lines)):
        lines[i] = lines[i][:-1]
        lines[i] = lines[i].split(' ')
        
        for j in range(0, len(lines[i])):
            lines[i][j] = float(lines[i][j])
            
        print(lines[i])

    return lines

def plot_cmap(leg, result):
    '''Helper function to plot temperature color maps based on specific points in time'''
    Na, Nb, Nc = len(result), len(result[0]), len(result[0][0])
    for i in range(Na):
        max_value, max_value_new = 0., 0.
        min_value, min_value_new = 0., 0.
        for j in range(Nb):
            max_value_new = max(result[i][j])
            
            if max_value_new > max_value:
                max_value = max_value_new
            elif min_value_new < min_value:
                min_value = min_value_new
            else: pass
        
        # Getting Jet colormap
        jet = cm.get_cmap('jet', 128)
        colormaps = [jet]
        
        cmap_len = len(colormaps)

        fig, axs = plt.subplots(1, cmap_len, figsize= ((cmap_len * 2) + 2, 3), constrained_layout= True, squeeze= False)
        plt.title(leg[i])
        for [ax, cmap] in zip(axs.flat, colormaps):
            psm = ax.pcolormesh(result[i], cmap= cmap, rasterized= True, vmin= max_value, vmax= min_value)
            fig.colorbar(psm, ax= ax)

def plot_quiver():
    plt.figure()
    plt.quiver(x,y,np.transpose(uplot),np.transpose(vplot)) 
    plt.title(f'Velocity quiver for time {t:.2f}')
    aux_arr = (u_star, v_star, p, u, v, psi, uplot, vplot)
    for i in range(len(result_params)):
        path = f'cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{reynolds[k]}/t_{t:.2f}/{i}_{result_params[i]}.txt'
        file = open(path, 'w')
        np.savetxt(path, np.transpose(aux_arr[i]))
        file.close()

def plot_values():
    values = [read_np_file(path + '/0_u_star.txt')]
    leg = ['u_star']
    plot_cmap(leg, values)
    
    values = [read_np_file(path + '/1_v_star.txt')]
    leg = ['v_star']
    plot_cmap(leg, values)
    
    values = [read_np_file(path + '/2_pressure.txt')]
    leg = ['pressure']
    plot_cmap(leg, values)
    
    values = [read_np_file(path + '/3_u.txt')]
    leg = ['u_new']
    plot_cmap(leg, values)
    
    values = [read_np_file(path + '/4_v.txt')]
    leg = ['v_new']
    plot_cmap(leg, values)
    
    values = [read_np_file(path + '/4_v.txt')]
    leg = ['v_new']
    plot_cmap(leg, values)

    plt.show()




if __name__ == "__main__":
    
    path = "/Users/felipeandrade/Documents/UnB/9. Semester/cavity-cfd/cavity_results/1.00x1.00"
    
    from os import walk

    f = []
    g = []
    h = []
    
    for (dirpath, dirnames, filenames) in walk(path):
        f.extend(filenames)
        g.extend(dirnames)
        h = dirpath.split('/')[1:]
        break
    
    print(f"{f}\n{g}\n{h}\n")
    
        
