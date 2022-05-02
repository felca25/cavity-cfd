import numpy as np
from matplotlib import pyplot as plt

se_fudeu = False


def calculate_u_star(u_star, u, v, Nx, Ny, dx, dy, dt, Re):
    if not se_fudeu:
        for j in range(0, Ny):
            for i in range(1, Nx):
                
                C1 = 0.25*(v[j+1,i] + v[j+1,i-1] + v[j,i] + v[j,i-1])
                
                R_conv = (-dt * ((u[j,i] * ((u[j,i+1] - u[j,i-1]) / (2*dx)))
                        +(C1 * ((u[j+1,i] - u[j-1,i]) / (2*dy)))))
                
                R_dif =  (dt/Re) * (((u[j,i+1] - 2*u[j,i] + u[j,i-1]) / (dx*dx))
                                +((u[j+1,i] - 2*u[j,i] + u[j-1,i]) / (dy*dy)))
                
                u_star[j,i] = u[j,i] + R_conv + R_dif
        
        
        for i in range(0, Nx+1):
            u_star[-1,i] = - u_star[0,i]
            u_star[Ny,i] = 2.0 - u_star[Ny-1,i]
            
        for j in range(0, Ny):
            u_star[j,0] = 0.0
            u_star[j,Nx] = 0.0

        return u_star

def calculate_v_star(v_star, u, v, Nx, Ny, dx, dy, dt, Re):
    if not se_fudeu:
        for j in range(1, Ny):
            for i in range(0, Nx):
                
                C2 = 0.25*(u[j,i+1] + u[j,i] + u[j-1,i+1] + u[j-1,i])
                
                R_conv = (-dt * ((v[j,i] * ((v[j,i+1] - v[j,i-1]) / (2*dx)))
                                +(C2 * ((v[j+1,i] - v[j-1,i]) / (2*dy)))
                                ))
                
                R_dif =  (dt/Re) * (((v[j,i+1] - 2*v[j,i] + v[j,i-1]) / (dx*dx))
                                +((v[j+1,i] - 2*v[j,i] + v[j-1,i]) / (dy*dy))
                                )
                
                v_star[j,i] = v[j,i] + R_conv + R_dif
                
        for j in range(0, Ny+1):
            v_star[j, -1] = -v_star[j,0]
            v_star[j,Nx] = -v_star[j, Nx-1]
            
        for i in range(0, Nx):
            v_star[0, i] = 0.0
            v_star[Ny, i] = 0.0

        return v_star

def calculate_pressure(p, u_star, v_star, Nx, Ny, dx, dy, dt, TOL=1e-5):
    global se_fudeu
    if not se_fudeu:
        erro = 10.
        iterations = 0
        while erro > TOL:
            R_max = 0
            for j in range(0, Ny):
                for i in range(0, Nx):
                    
                    if i == 0 and j == 0:
                        alpha = -((1. / (dx*dx)) + (1. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = ((p[j, i+1] - 1.*p[j,i]) / (dx*dx))
                        R2_y = ((p[j+1, i] - 1.*p[j,i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y
                            
                    elif i == 0 and j == Ny-1:
                        
                        alpha = -((1. / (dx*dx)) + (1. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = ((p[j, i+1] - 1.*p[j,i]) / (dx*dx))
                        R2_y = ((- 1.*p[j,i] + p[j-1, i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y
                        
                    elif i == 0 and j != 0 and j != Ny-1:
                                        
                        alpha = -((1. / (dx*dx)) + (2. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = ((p[j, i+1] - 1.*p[j,i]) / (dx*dx))
                        R2_y = ((p[j+1, i] - 2.*p[j,i] + p[j-1, i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y
                    
                    elif i == Nx-1 and j == 0:
                            
                        alpha = -((1. / (dx*dx)) + (1. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = -((- 1.*p[j,i] + p[j, i-1]) / (dx*dx))
                        R2_y = -((p[j+1, i] - 1.*p[j,i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y
                            
                    elif i == Nx-1 and j == Ny-1:
                            
                        alpha = -((1. / (dx*dx)) + (1. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = -((- 1.*p[j,i] + p[j, i-1]) / (dx*dx))
                        R2_y = -((- 1.*p[j,i] + p[j-1, i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y
                        
                    elif i == Nx-1 and j != 0 and j != Ny-1:
                        alpha = -((1. / (dx*dx)) + (2. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = ((p[j, i-1] - 1.*p[j,i]) / (dx*dx))
                        R2_y =-((p[j+1, i] - 2*p[j,i] + p[j-1, i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y
                    
                    elif i != 0 and i != Nx-1 and j == 0:
                            
                        alpha = -((2. / (dx*dx)) + (1. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = ((p[j, i+1] - 2*p[j,i] + p[j, i-1]) / (dx*dx))
                        R2_y = ((p[j+1, i] - 1.*p[j,i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y
                            
                    elif i != 0 and i != Nx-1 and j == Ny-1:
                            
                        alpha = -((2. / (dx*dx)) + (1. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = ((p[j, i+1] - 2*p[j,i] + p[j, i-1]) / (dx*dx))
                        R2_y = ((p[j-1, i] - 1.*p[j,i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y

                    else:
                        alpha = -((2. / (dx*dx)) + (2. / (dy*dy)))
                        
                        R1_x = ((u_star[j, i+1] - u_star[j, i]) / (dt*dx))
                        R1_y = ((v_star[j+1, i] - v_star[j, i]) / (dt*dy))
                        
                        R2_x = ((p[j, i+1] - 2*p[j,i] + p[j, i-1]) / (dx*dx))
                        R2_y = ((p[j+1, i] - 2*p[j,i] + p[j-1, i]) / (dy*dy))
                
                        Residual = R1_x + R1_y - R2_x - R2_y
                        
                    Residual /= alpha 
                    p[j,i] = p[j,i] + (Residual / alpha)
                    
                    if np.abs(p[j,i]) >= 1000:
                        se_fudeu = True
                        print(f'Element [{i}][{j}]')
                        print(f"R = {Residual}, {R1_x}, {R1_y}, {R2_x}, {R2_y}")
                        raise ValueError
                
                    if np.abs(Residual) >= R_max:
                        R_max = np.abs(Residual)
                    
            erro = R_max
            iterations += 1
            
            # print(f"Iteration {iterations}, erro = {erro}")
        
        for i in range(0,Nx):
            p[-1,i] = p[0,i]
            p[Ny,i] = p[Ny-1,0]
            
        for j in range(0, Ny):
            p[j,-1] = p[j,0]
            p[j,Nx] = p[j, Nx-1]
            
        p[-1,-1] = p[0,0]
        p[Ny,-1] = p[Ny-1,0]
        p[-1,Nx] = p[0,Nx-1]
        p[Ny,Nx] = p[Ny-1,Nx-1]
        
        # print(iterations)
        return p

def calculate_u_new(u, u_star, p, Nx, Ny, dx):
    if not se_fudeu:
        for j in range(-1, Ny+1): 
            for i in range(1, Nx):
                u[j,i] = u_star[j,i] - (dt*((p[j,i] - p[j,i-1]) / dx))
        
        return u

def calculate_v_new(v, v_star, p, Nx, Ny, dx):
    if not se_fudeu:
        for j in range(1, Ny):
            for i in range(-1, Nx+1):
                v[j,i] = v_star[j,i] - (dt*((p[j,i] - p[j-1,i]) / dy))
        
        return v
                    


if __name__ == "__main__":
    Nx = 5
    Ny = 5
    
    Re = 10
    
    dx = 1 / Nx
    dy = 1 / Ny
    
    dt = 0.25 * Re * (dx**2)
    
    if dx >= 1/np.sqrt(Re):
        print('Se fudeu')
    
    TOL = 1e-6
    
    u = np.zeros([Ny+2, Nx+1], float)
    v = np.zeros([Ny+1, Nx+2], float)
    p = np.zeros([Ny+2, Nx+2], float)
    
    u[Ny, :] = 2.0
    
    u_star = np.copy(u)
    v_star = np.copy(v)
    
    time = 0.0
    k = 0
    while k <= 200:
        print(k)
        if not se_fudeu:
            u_star = calculate_u_star(u_star, u, v, Nx, Ny, dx, dy, dt, Re)
            v_star = calculate_v_star(v_star, v, u, Nx, Ny, dx, dy, dt, Re)
            p = calculate_pressure(p, u_star, v_star, Nx, Ny, dx, dy, dt, TOL)
            u = calculate_u_new(u, u_star, p, Nx, Ny, dx)
            v = calculate_v_new(v, v_star, p, Nx, Ny, dx)
            k +=1
        else:
            break
    
    print(u_star[::-1])
    print(v_star[::-1])
    print(p[::-1])
    print(u[::-1])
    print(v[::-1])
    
