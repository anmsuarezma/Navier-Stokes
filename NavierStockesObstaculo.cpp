//Fludo estacionario- viscoso - incompresible 
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>

const int nx = 101; // Pueden ir hasta 1001 con la velocidad 10.0. Si aumenta la velocidad
const int ny = 101; //también puede aumentar el a nx y ny
const double lx = 10.;
const double ly = 10.;
const double dx = lx/(nx-1);
const double dy = ly/(ny-1);
const double V0 = 5.;
const double R = 2.;
const double nu = dx*V0/R;
const double rel = 0.1;

class NavierStockes
{

private:

    double phi[nx+2][ny+2];
    double vort[nx+2][ny+2];

public:
    void Flujo(void);
    void condicionesDeFronteraFlujo(void);
    void imprimase(void);
    double fuenteExterna(void);
    void condicionesIniciales(void);
    void gnuplot(void);
};

void NavierStockes::condicionesDeFronteraFlujo(void)
{
    double dx2=dx*dx;
    double dx2Inv=1./dx2;
    double a1 = 1.-rel;   
    double newPhi[nx+2][ny+2] = {};
    double newVort[nx+2][ny+2] = {};
    #pragma omp parallel for schedule(auto)
    for (int i=0; i<nx; ++i)
        for (int j=0; j<ny; ++j)
            newPhi[i][j] = phi[i][j];

    #pragma omp parallel for schedule(auto)
    for (int i=0; i<nx; ++i)
        for (int j=0; j<ny; ++j)
            newVort[i][j] = vort[i][j];

    //Abajo
    #pragma omp parallel for 
    for (int i=1; i<nx-1; i++)
    {
        phi[i][0] = 0.0;   
        vort[i][0] = 0.0;     
    }

    //#pragma omp parallel for 
    //for (int i=1; i<nx-1; i++)
      //  vort[i][0] = 0.0;     

    //Arriba
    #pragma omp parallel for 
    for (int i=1 ; i<nx-1; i++)
    { 
        vort[i][ny-1] = 0.0;
        phi[i][ny-1] = a1*newPhi[i][ny-1] + 0.25*rel*(2.*phi[i][ny-2] + 2.*V0*dy + phi[i-1][ny-11] + phi[i+1][ny-1] - dx2*vort[i][ny-1]);  
    }
   // #pragma omp parallel for 
    //for (int i=1; i<nx-1; i++)
    	//phi[i][ny-1] = a1*newPhi[i][ny-1] + 0.25*rel*(2.*phi[i][ny-2] + 2.*V0*dy + phi[i-1][ny-11] + phi[i+1][ny-1] - dx2*vort[i][ny-1]);  

    //Izquierda
    #pragma omp parallel for 
    for (int j=2; j<ny-2; j++)
    {
        phi[0][j] = a1*newPhi[0][j] + 0.25*rel*(2.*phi[1][j] + phi[0][j+1] + phi[0][j-1] - dx2*vort[0][j]);  
        vort[0][j] = 0.0;     
    }
        
    //#pragma omp parallel for 
    //for (int j=2; j<ny-2; j++)
        //vort[0][j] = 0.0;      

    //Derecha
    #pragma omp parallel for 
    for (int j=1; j<ny-1; j++)
    {
        phi[nx-1][j] = a1*newPhi[nx-1][j] + 0.25*rel*(2.*phi[nx-2][j] + phi[nx-1][j-1] + phi[nx-1][j+1] - dx2*vort[nx-1][j]); 
        vort[nx-1][j] = a1*newVort[nx-1][j] + 0.25*rel*(2.*vort[nx-2][j] + vort[nx-1][j-1] + vort[nx-1][j+1]);
    }
        
        
    //#pragma omp parallel for 
    //for (int j=2; j<ny-2; j++)
        //vort[nx-1][j] = a1*newVort[nx-1][j] + 0.25*rel*(2.*vort[nx-2][j] + vort[nx-1][j-1] + vort[nx-1][j+1]);

    //Obstaculo
    #pragma omp parallel for 
    for (int j=1; j<ny/3; ++j)
    {
        phi[nx/3][j] = 0.; 
        phi[nx/3 + nx/4][j] = 0.;
        vort[nx/3][j] = (2.*dx2Inv)*(phi[nx/3-1][j] - phi[nx/3][j]); //pared izquierda
    }
        
   /* #pragma omp parallel for         
    for (int j=1; j<ny/3; ++j)
        phi[nx/3 + nx/4][j] = 0.;*/
            
    #pragma omp parallel for 
    for (int i=nx/3; i<nx/3 + nx/4; ++i)
        phi[i][ny/3] = 0.;

   /* #pragma omp parallel for 
    for (int j=1; j<ny/3; ++j)
        vort[nx/3][j] = (2.*dx2Inv)*(phi[nx/3-1][j] - phi[nx/3][j]); //pared izquierda
*/
    #pragma omp parallel for 
    for (int j=1; j<ny/3; ++j)
        vort[nx/3+ nx/4][j] = (2.*dx2Inv)*(phi[nx/3+ nx/4 + 1][j] - phi[nx/3 + nx/4][j]); //pared derecha
    
    /*#pragma omp parallel for 
    for (int i=nx/3; i<nx/3 + nx/4; ++i)
        vort[i][ny/3] = (2.*dx2Inv)*(phi[i][ny/3 + 1] - phi[i][ny/3]); //pared superior
    */
   #pragma omp parallel for 
    for (int ii=nx/3; ii<nx/3 + nx/4 ; ++ii)
    {
        vort[ii][ny/3] = (2.*dx2Inv)*(phi[ii][ny/3 + 1] - phi[ii][ny/3]); //pared superior
        for (int jj=0; jj<ny/3-1; ++jj)
        {
            vort[ii][jj] = 0.;
            phi[ii][jj] = 0.;
        }
    }
    
   /*#pragma omp parallel for 
    for (int ii=nx/3; ii<nx/3 + nx/4 ; ++ii)
        for (int jj=0; jj<ny/3-1; ++jj)
            phi[ii][jj] = 0.;*/
}


void NavierStockes::condicionesIniciales(void)
{
        
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; ++j)
            phi[i][j] = 1.*j*dy;

    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; ++j) 
            vort[i][j] = 0.0;
}

void NavierStockes::Flujo(void)
{
    double dx2=dx*dx;
    double a1 = 1.-rel;   
    double a2 = 0.25/nu;
    double newPhi[nx][ny] = {};
    double newVort[nx][ny] = {};
    double Aux1[nx][ny] = {};
    double Aux2[nx][ny] = {};
    double Aux3[nx][ny] = {};

    #pragma omp parallel for
    for (int i=0; i<nx; ++i)
    {
        for (int j=0; j<ny; ++j)
        {
            newPhi[i][j] = phi[i][j];
            newVort[i][j] = vort[i][j]; 
        }        
    }

   /* #pragma omp parallel for
    for (int i=0; i<nx; ++i)
        for (int j=0; j<ny; ++j)
            newVort[i][j] = vort[i][j];
    */         
      
    #pragma omp parallel for
    for (int i=1; i<nx-1; ++i)
        for (int j=1; j<ny-1; ++j)
            phi[i][j] = a1*newPhi[i][j] + 0.25*rel*(newPhi[i][j+1] + newPhi[i][j-1] + newPhi[i+1][j] + newPhi[i-1][j]- dx2*newVort[i][j]);
            //phi[i][j] = (0.25*rel)*(newPhi[i+1][j] + newPhi[i-1][j] + newPhi[i][j+1] + newPhi[i][j-1] - vort[i][j]) + (1. - rel)*newPhi[i][j];

    #pragma omp parallel for
    for (int i=1; i<nx-1; ++i)
    {
        for (int j=1; j<ny-1; ++j)
        {
            Aux1[i][j] = newVort[i][j+1] + newVort[i][j-1] + newVort[i+1][j] + newVort[i-1][j];
            Aux2[i][j] = (phi[i][j+1] - phi[i][j-1])*(newVort[i+1][j] - newVort[i-1][j]);
            vort[i][j] = a1*newVort[i][j] + 0.25*rel*(Aux1[i][j] + a2*(-Aux2[i][j] + Aux3[i][j]));
        }       
    }

    /*#pragma omp parallel for
    for (int i=1; i<nx-1; ++i)
        for (int j=1; j<ny-1; ++j)
            Aux2[i][j] = (phi[i][j+1] - phi[i][j-1])*(newVort[i+1][j] - newVort[i-1][j]);


    #pragma omp parallel for
    for (int i=1; i<nx-1; ++i)
        for (int j=1; j<ny-1; ++j)        
            Aux3[i][j] = (phi[i+1][j] - phi[i-1][j])*(newVort[i][j+1] - newVort[i][j-1]);
    */
     
    #pragma omp parallel for
    for (int i = 1; i < nx-1; ++i)
        for (int j = 1; j < ny-1; ++j)
            vort[i][j] = a1*newVort[i][j] + 0.25*rel*(Aux1[i][j] + a2*(-Aux2[i][j] + Aux3[i][j]));

}

double NavierStockes::fuenteExterna(void)
{
    double F = 1.0;
    return F;    
}

void NavierStockes::imprimase(void)
{
     std::cout << "splot '-' w l lw 2 " << std::endl;

    double x,y;

    for (int ii=0; ii<nx; ++ii)
    {
        x = ii * dx;
        for (int jj=0; jj<ny; ++jj)
        {
            y = jj*dy;
            std::cout << x << " \t" << y << " \t" << phi[ii][jj] << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "e" << std::endl;
}


void NavierStockes::gnuplot(void)
{
    std::cout << "set terminal gif animate delay 100" << std::endl;
    std::cout << "set out 'NSObstaculos.gif'" << std::endl;
    //std::cout << "set contour base" << std::endl;
    std::cout << "set pm3d" << std::endl;
    std::cout << "unset surface" << std::endl; // Comentar para ver 3D
    //std::cout << "set cntrparam levels auto 10" << std::endl;
    std::cout << "set view map" << std::endl; // Comentar para ver 3D
    std::cout << "set palette rgb 21,22,23" << std::endl;// Comentar para ver 3D
    std::cout << "set xlabel 'x' " << std::endl;
    std::cout << "set ylabel 'y' " << std::endl;
    std::cout << "unset key" << std::endl;

}

int main(void)
{
    NavierStockes fluido;
    int nit = 1000;
    fluido.condicionesIniciales();
    fluido.condicionesDeFronteraFlujo();
    fluido.gnuplot();
    //#pragma omp parallel for schedule(auto)
    for (int t = 0; t <= nit; t++)
    {
        fluido.Flujo();
        fluido.condicionesDeFronteraFlujo();
        fluido.imprimase();
    }  

    //fluido.imprimase();
    return 0;
}
