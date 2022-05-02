#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

const int Nx = 5;
const int Ny = 5;

const float Lx = 1.0;
const float Ly = 1.0;

const int Re = 10;

const float dx = Lx / Nx;
const float dy = Ly / Ny;

const float dt = 0.25 * Re * (dx * dx);

const float TOL = 1e-6;

float u[Nx+1][Ny+2] = {};
float v[Nx+1][Ny+2] = {};
float p[Nx+2][Ny+2] = {};

float u_star[Nx+1][Ny+2] = {};
float v_star[Nx+2][Ny+1] = {};

int main(){
    for(int i = 0; i < Nx; i++){
        u[i][Ny] = 2.0;
    }


    // Calculating u_star
    for(int i = 1; i < Nx; i++){
        for(int j = 0; j < Nx; j++){
            float C1 = 0.25 * (v[i][j+1] + v[i-1][j+1] + v[i][j] + v[i-1][j]);
            float R_conv_x = - dt * ((u[i][j] * (u[i+1][j] - u[i-1][j])) / (2*dx));
            float R_conv_y = - dt * ((C1 * (u[i+1][j] - u[i-1][j])) / (2*dy));
            float R_dif = (dt / Re) * ((u[i+1][j] - 2*u[i][j] + u[i-1][j]) / (dx*dx))
                                    + ((u[i][j+1] - 2*u[i][j] + u[i][j-1]) / (dy*dy));
            u_star[i][j] = u[i][j] + R_conv_x + R_conv_y + R_dif;
        }
    }
        // Boundary conditions
    for(int i = 0; i < Nx+1; i++){
        u_star[i][Ny+1] = - u_star[i][0];
        u_star[i][Ny] = 2.0 - u_star[i][Ny-1];
    }

    for(int j = 0; j < Ny; j++){
        u_star[0][j] = 0.0;
        u_star[Nx][j] = 0.0;
    }

    // Calculating v_star
    for(int i = 0; i < Nx; i++){
        for(int j = 1; j < Nx; j++){
            float C2 = 0.25 * (u[i+1][j] + u[i][j] + u[i+1][j-1] + u[i][j-1]);
            float R_conv_x = - dt * ((C2 * (v[i+1][j] - v[i-1][j])) / (2*dx));
            float R_conv_y = - dt * (v[i+1][j] * (v[i+1][j] - v[i-1][j])) / (2*dy);
            float R_dif = (dt / Re) * ((v[i+1][j] - 2*v[i][j] + v[i-1][j]) / (dx*dx))
                                    + ((v[i][j+1] - 2*v[i][j] + v[i][j-1]) / (dy*dy));
            v_star[i][j] = v[i][j] + R_conv_x + R_conv_y + R_dif;
        }
    }
        // Boundary conditions
    for(int j = 0; j < Ny+1; j++){
        v_star[Nx+1][j] = - v_star[0][j];
        v_star[Nx][j] = - v_star[Ny-1][j];
    }

    for(int i = 0; i < Nx; i++){
        u_star[i][0] = 0.0;
        u_star[i][Ny] = 0.0;
    }

    // Calculating pressure
    float error = 10;
    float alpha;
    float R1_x;
    float R1_y;
    float R2_x;
    float R2_y;
    float Residual;

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j ++){
            if (i == 0 && j == 0){
                alpha = - ((1.0/(dx*dx) + (1.0/(dy*dy))));
                R1_x = ((u_star[i+1][j] - u_star[i][j]) / (dt*dx));
                R1_y = ((v_star[i][j+1] - v_star[i][j]) / (dt*dy));

                R2_x = ((p[i+1][j] - p[i][j]) / (dx*dx));
                R2_y = ((p[i][j+1] - p[i][j]) / (dy*dy));
                Residual = R1_x + R1_y - R2_x - R2_y;

            } else if (i == 0 && j == Ny-1){
                alpha = - ((1.0/(dx*dx) + (1.0/(dy*dy))));
                R1_x = ((u_star[i+1][j] - u_star[i][j]) / (dt*dx));
                R1_y = ((v_star[i][j+1] - v_star[i][j]) / (dt*dy));

                R2_x = ((p[i+1][j] - p[i][j]) / (dx*dx));
                R2_y = ((-p[i][j] + p[i][j-1]) / (dx*dy));
                Residual = R1_x + R1_y - R2_x - R2_y;
            } else if (i == Nx-1 && j == 0){
                
            }
        }
    }







    for(int i = 0; i < Nx+1; i++){
        for(int j = 0; j < Ny+2; j++){
            cout << u_star[i][j] << ' ';
        }
        cout << "\n";
    }
    for(int i = 0; i < Nx+2; i++){
        for(int j = 0; j < Ny+1; j++){
            cout << v_star[i][j] << ' ';
        }
        cout << "\n";
    }
    for(int i = 0; i < Nx+2; i++){
        for(int j = 0; j < Ny+2; j++){
            cout << v_star[i][j] << ' ';
        }
        cout << "\n";
    }


    return 0;
}