class Cavity:
    
    def __init__(self, Nx, Ny, Lx,Ly, Re, t_final, TOL=1e-5):
        
        self.Nx, self.Ny = Nx, Ny
        self.Lx, self.Ly = Lx, Ly
        self.Re = Re
        self.dx, self.dy = self.Lx / self.Nx, self.Ly / self.Ny
        self.dt = self.reset_dt()
        self.t_final = t_final
        self.TOl = TOL
        
        self.create_paths()/


    def reset_dt(self):
        
        self.dt = 0.25 * self.Re * (self.dx**2)
        
        while self.dx >= 1 / (self.Re ** (1/2)):
            self.dx = (1 / (self.Re ** (1/2))) - 1e-3
        
        while self.dt >= self.dx:
            self.dt = self.dx - 1e-3

        self.t_final = 0.5
        
        return self.dt
    
    def create_paths(Re, t_arr, result_params):
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
    
    for i in range(len(Re)):
        for j in range(len(t_arr)):
            
            path = f'cavity_results/{Lx:.2f}x{Ly:.2f}/Re_{Re[i]}/t_{t_arr[j]:.2f}'
            
            try:
                os.makedirs(path)
                
            except FileExistsError:
                print(f'{path} Already Exists')
                
            for k in range(len(result_params)):
                txt_path = f'{path}/{k}_{result_params[k]}.txt'
                f =  open(txt_path, 'w').close()
                paths.append(txt_path)
                
    return paths
