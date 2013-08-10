
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//==================================================
// structure to hold parameters for model
//==================================================
typedef struct {
    double dt;
    double anisotropy;
    double latent_heat;
    double w_phi;
} Parameters;

//==================================================
double ** alloc_2d(int nx, int ny)

    // return 2d array allocated with nx by ny elements
//==================================================
{
    int i;
    double ** array = (double **) malloc(sizeof(double *)*nx);
    for (i=0; i<nx; i++)
        array[i] = (double *) malloc(sizeof(double)*ny);
    return array;
}

//==================================================
void free_2d(double ** array, int nx)
    // free 2d array
//==================================================
{
    int i;
    for (i=0; i<nx; i++)
        free(array[i]);
    free(array);
}

//==================================================
void init_droplet(double ** phi, double ** temp, int nx, int ny)

    // initialize small circular droplet in center of the system
    // droplet is at the melting temperature
    // the surrounding material is undercooled
//==================================================
{
    int i,j;
    double droplet_radius = 10;

    for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
        phi[i][j] = -1;
        temp[i][j] = 0.5; // temperature of surrounding material
    }

    for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
    {
        double rx = i - nx/2.0;
        double ry = j - ny/2.0;
        double r = sqrt(rx*rx + ry*ry);

        if (r < droplet_radius) {
            temp[i][j] = 1;
            phi[i][j] = 1;
        }
    }
}

//==================================================
void output_tecplot(const char * filename,  /* in */
                    const char * mode,      /* in */
                    double ** phi,          /* in */
                    double ** temp,         /* in */
                    int nx, int ny,         /* in */
                    double time)            /* in */

    /* Output order parameter and temperature field 
    in tecplot animation format */
//==================================================
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

//==================================================
inline double laplacian_5pt(double **data,  /* in */
                            int i, int j)   /* in */

    // returns numerator of discretized 5-point
    // laplacian operator about the gridpoint (i,j)
//==================================================
{
    double nabla_sq = data[i-1][j] +
                      data[i+1][j] +
                      data[i][j-1] +
                      data[i][j+1] -
                      4*data[i][j];
    
    return nabla_sq;
}

//==================================================
inline double div_w_grad_phi(double ** w,    /* in */
                             double ** phi,  /* in */
                             int i, int j)   /* in */

    // returns divergence of the product of the 
    // first array with the gradient of the second
    // at location (i,j) in grid
//==================================================
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

//==================================================
inline double grad_x( double ** array,   /* in */
                      int i, int j )     /* in */

    // returns centered difference in x-dimensions
//==================================================
{
    return 0.5*(array[i+1][j] - array[i-1][j]);
}

//==================================================
inline double grad_y( double ** array,   /* in */
                      int i, int j )     /* in */

    // returns centered difference in y-dimensions
//==================================================
{
    return 0.5*(array[i][j+1] - array[i][j-1]);
}

//==================================================
inline double dg(double p)
    // returns the derivative of the double-well potential
//==================================================
{
    return p*p*p - p;
}


//==================================================
inline double dp(double p)
    // returns derivative of the solidification potential offset
//==================================================
{
    double q = (1.0 - p*p);
    return q*q;
}

//==================================================
void solve(double ** phi,        /* in-out */
           double ** temp,       /* in-out */
           int nx, int ny,       /* in */
           Parameters * params)  /* in */

    // Integrate evolution equation for order parameter
    // and the heat equation for the temperature field.
    // All required parameters contained within parameter struct.
    // Two ghost rows required on all sides of system
//==================================================
{

    int i,j;
    double tau[nx][ny];
    double w[nx][ny];
    double w_prime[nx][ny];
    double ** w_sq = alloc_2d(nx, ny);

    // loop thru system, including the inside ghost row
    for (i=1; i<nx-1; i++)
    for (j=1; j<ny-1; j++)
    {
        // calculate anisotropic surface energy parameter

        double phi_x = grad_x(phi, i, j);
        double phi_y = grad_y(phi, i, j);

        double theta = atan(phi_y/(phi_x+0.000001));
        double A = 1.0 + params->anisotropy * cos(4.0*theta);

        w[i][j] = params->w_phi * A;
        w_prime[i][j] = - params->w_phi * 4.0 * params->anisotropy * sin(4.0*theta);
        w_sq[i][j] = w[i][j]*w[i][j];
        tau[i][j] = A*A;
    }

    // loop thru system again, not including ghost rows
    for (i=2; i<nx-2; i++)
    for (j=2; j<ny-2; j++)
    {
        double p = phi[i][j];

        // use properly centered derivatives for terms from chain rule
        // p850 Numerical Recipes C
        double term1 = w[i+1][j] * w_prime[i+1][j] * grad_y(phi, i+1, j) -
                       w[i-1][j] * w_prime[i-1][j] * grad_y(phi, i-1, j);

        double term2 = w[i][j+1] * w_prime[i][j+1] * grad_x(phi, i, j+1) -
                       w[i][j-1] * w_prime[i][j-1] * grad_x(phi, i, j-1);

        double interface_energy = div_w_grad_phi(w_sq, phi, i, j);

        // chemical potential calculated straight from equation in the book
        double dphi_dt = interface_energy - term1 + term2 
                         - dg(p) - params->latent_heat*(temp[i][j] - 1.0)*dp(p);

        // update order parameter
        phi[i][j] += params->dt * dphi_dt / tau[i][j];

        // update temperature
        temp[i][j] += params->dt * ( laplacian_5pt(temp, i, j) + 0.5*dphi_dt );
    }

    free_2d(w_sq, nx);
}

int main()
{

    int nx = 400;
    int ny = 400;
    int nsteps =  10000;
    int out_freq = 1000;

    Parameters * params = (Parameters *) malloc(sizeof(Parameters));
    params->dt          = 0.2;
    params->anisotropy  = 0.05;
    params->latent_heat = 3.0;
    params->w_phi       = 1.0;

    int istep;
    double ** phi = alloc_2d(nx, ny);
    double ** temp = alloc_2d(nx, ny);

    init_droplet(phi, temp, nx, ny);
    output_tecplot("out.dat", "w", phi, temp, nx, ny, 0);

    for (istep=0; istep<nsteps; istep++)
    {
        // iterate evolution equations
        solve(phi, temp, nx, ny, params);

        // no need to do boundary conditions
        // unless dendrites reach each of system

        // write data to output file
        if ((istep+1)%out_freq==0) 
        {
            double percent = (istep+1)/((double)nsteps)*100;
            printf("%d: %5.0f%% \n", istep+1, percent);

            output_tecplot("out.dat", "a", phi, temp, nx, ny, istep+1);
        }
    }

    free_2d(phi, nx);
    free_2d(temp, nx);

    return 0;

}

