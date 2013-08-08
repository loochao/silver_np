
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "barr.h"

void init_droplet(double ** phi, double ** temp, int nx, int ny)
{
    int i,j;

    for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
        phi[i][j] = -1;
        temp[i][j] = 0.5;
    }

    for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
        double rx = i - nx/2;
        double ry = j - ny/2;
        double r = sqrt(rx*rx + ry*ry);

        if (r<10) {
            temp[i][j] = 1;
            phi[i][j] = 1;
        }
    }
    for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
        double rx = i - nx/2+10;
        double ry = j - ny/2+10;
        double r = sqrt(rx*rx + ry*ry);

        if (r<10) {
            temp[i][j] = 1;
            phi[i][j] = 1;
        }
    }
}

/* Output tecplot format */
void output_tecplot(const char * filename, const char * mode, double ** phi, double ** temp, int nx, int ny, double time)
{
    int i,j;
    FILE * fp = fopen(filename, mode);

    fprintf(fp, "VARIABLES = X Y phi temp\n");
    fprintf(fp, "ZONE I = %d, J = %d, strandid=1, solutiontime=%f\n", nx, ny, time);

    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
        fprintf(fp, "%d %d %f %f\n", i, j, phi[i][j], temp[i][j]);

    fclose(fp);
}

inline double laplacian_5pt(double ** data, int i, int j)
{
    double nabla_sq = data[i-1][j] +
                      data[i+1][j] +
                      data[i][j-1] +
                      data[i][j+1] -
                      4.0*data[i][j];

    return nabla_sq;
}

inline double grad_sq(double ** data, int i, int j)
{
    double grad_x = (data[i+1][j] - data[i-1][j]);
    double grad_y = (data[i][j+1] - data[i][j-1]);
    return (grad_x*grad_x + grad_y*grad_y)*0.25;
}

inline double div_w_grad_phi(double ** w,
                             double ** phi,
                             int i, int j)
{

        double mx1 = w[i+1][j] + w[i][j];
        double mx2 = w[i][j] + w[i-1][j];

        double my1 = w[i][j+1] + w[i][j];
        double my2 = w[i][j] + w[i][j-1];

        double term1 = mx1 * (phi[i+1][j] - phi[i][j]) -
                       mx2 * (phi[i][j] - phi[i-1][j]);

        double term2 = my1 * (phi[i][j+1] - phi[i][j]) -
                       my2 * (phi[i][j] - phi[i][j-1]);


        double ddt = 0.5 * (term1 + term2);

        return ddt;

}    

inline double dg(double p)
{
    return p*p*p - p;
}


inline double dp(double p)
{
    double q = (1.0 - p*p);
    return q*q;
}

void solve(double ** phi, double ** temp, int nx, int ny)
{

    int i,j;
    double dt = 0.01;
    double eps = 0.5;
    int dims[2] = {nx, ny};

    double chem_potential[nx][ny];
    double lap_temp[nx][ny];
    double ** w_sq = barr_alloc(sizeof(double), 2, dims, 0);
    double tau[nx][ny];

    for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
        double p = phi[i][j];
        double nabla_sq = laplacian_5pt(phi, i, j);    
        double r = ((double)rand())/((double)RAND_MAX)*2-1;
        //chem_potential[i][j] = dg(p) + (temp[i][j]-1)*dp(p) - nabla_sq + 0.5*r;
        lap_temp[i][j] = laplacian_5pt(temp, i, j);

        double grad_x = (phi[i+1][j] - phi[i-1][j]) * 0.5;
        double grad_y = (phi[i][j+1] - phi[i][j-1]) * 0.5;

        double theta = atan(grad_y/(grad_x+0.000001));
        double A = 1.0 + eps * cos(4.0 * theta);

        printf("%d %d %f %f %f %f %f %f\n", i, j, phi[i+1][j], phi[i-1][j], grad_x, grad_y, theta, A);
        w_sq[i][j] = A*A;
        tau[i][j] = A*A;
    }

    for (i=1; i<nx-1; i++)
    for (j=1; j<ny-1; j++)
    {
        double p = phi[i][j];
        phi[i][j] += dt * ( div_w_grad_phi(w_sq, phi, i, j) - dg(p) - (temp[i][j]-1)*dp(p));
        temp[i][j] += dt * (lap_temp[i][j] -chem_potential[i][j]);
    }

}

void free_energy(double ** phi, int nx, int ny, double *total, double * local, double * interfacial)
{
    int i,j;
    *total = 0;
    *local = 0;
    *interfacial = 0;
    for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
        double p = phi[i][j];
        double p_sq = p*p;
        double p_fr = p_sq*p_sq;

        *interfacial += 0.5 * grad_sq(phi,i,j);
        *local += -0.5*p_sq + 0.25*p_fr + p-2*p*p_sq/3.0 + p_fr*p/5.0;
        *total += *local + *interfacial;
    }
}

int main()
{

    int nx = 100;
    int ny = 100;
    int dims[2] = {nx, ny};

    int nsteps =  0;
    int out_freq = 1;
    int istep;
    double ** phi = (double **) barr_alloc(sizeof(double), 2, dims, 2);
    double ** temp = (double **) barr_alloc(sizeof(double), 2, dims, 2);

    init_droplet(phi, temp, nx, ny);
    barr_pbc(&(phi[0][0]), 2, dims, 2);
    barr_pbc(&(temp[0][0]), 2, dims, 2);
    output_tecplot("out.dat", "w", phi, temp, nx, ny, 0);

    for (istep=0; istep<nsteps; istep++)
    {

        solve(phi, temp, nx, ny);
        barr_pbc(&(phi[0][0]), 2, dims, 2);
        barr_pbc(&(temp[0][0]), 2, dims, 2);

        if ((istep+1)%out_freq==0) 
        {
            /*
            double total,local,interfacial;
            free_energy(phi, nx, ny, &total, &local, &interfacial);
            printf("%d %f %f %f\n", istep+1, total, local, interfacial);
            */
            printf("%d\n", istep+1);
            output_tecplot("out.dat", "a", phi, temp, nx, ny, istep+1);
        }
    }

    barr_free((double*)phi, 2);

    return 0;

}

