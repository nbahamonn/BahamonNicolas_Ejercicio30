#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// Codigo Modificado del Ejercicio 29.

// variable globales
const double L = 1.0;
const double rho = 0.01;
const int ten = 40;
const double Xi = 0.0;
const double Xf = 1.0;
const double Ti = 0.0;
const double Tf = 0.1;
const double dx = 0.01;
const double dt = 0.0001;
const int Nt = Tf/dt;
const int Nx = Xf/dx;
const double c = sqrt(ten/rho);
const double cl = dx/dt;
const double ratio = c*c/(cl*cl);



//Inicializacion 
int main(void){
    
    double **G = new double *[Nt+1];
        for(int i = 0; i <= Nt; i++){
            G[i] = new double[Nx+1];
        }
    
    for(int j = 0; j <= Nx; j++){
		G[0][j] = (10E-4)*sin((2*M_PI*j*dx)/L);
	}
    
    for(int j = 1; j < Nx; j++){
		G[1][j] = G[0][j] + ((c*c)/(2*cl*cl))*(G[0][j+1] + G[0][j-1] - 2*G[0][j]);
    }
    
    for(int i = 1; i < Nt; i++){
		for(int j = 1; j < Nx; j++){
				G[i+1][j] = 2*G[i][j] - G[i-1][j]+ ((c*c)/(cl*cl))*(G[i][j+1] + G[i][j-1] - 2*G[i][j]);
		}
	}
    
    ofstream outfile;
    outfile.open("Datos.dat");
	
	for(int i = 0; i <= Nt; i++){
		for(int j = 0; j <= Nx; j++){
			outfile << G[i][j] << "\t";
		}
		outfile << endl;
	}		
	outfile.close();
    return 0;
}