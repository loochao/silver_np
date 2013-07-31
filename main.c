#include <stdlib.h>
#include <stdio.h>
#include "solve.h"

double ** alloc_2d(int nx, int ny)
{
    int i;
    double ** array = malloc(sizeof(double) * nx);
    for (i=0; i<nx; i++)
        array[i] = malloc(sizeof(double) * ny);
    return array;
}

void free_2d(double ** array, int nx)
{
    int i;
    for (i=0; i<nx; i++)
        free(array[i]);
    free(array);
}

/* Output tecplot format */
void output_tecplot(const char * filename, const char * mode, double ** phi, int nx, int ny, double time)
{
    int i,j;
    FILE * fp = fopen(filename, mode);

    fprintf(fp, "VARIABLES = phi\n");
    fprintf(fp, "ZONE I = %d, J = %d, strandid=1, solutiontime=%f\n", nx, ny, time);

    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
        fprintf(fp, "%d %d %f\n", i, j, phi[i][j]);

    fclose(fp);
}

/* Set initial values of the phi array */
void initial_conditions(double ** phi, int nx, int ny)
{
  int i, j;
  for (i=0; i<nx; i++)
    {
      for(j=0; j<ny; j++)
        phi[i][j] = 0;
    }
  phi[nx/2][ny/2] = 1;
}

/* Apply Periodic Boundary Conditions */
void pbc(double ** phi, int nx, int ny)
{
  int i, j;
  for (i=0; i<nx; i++)
    phi[i][0] = phi[i][ny];

  for (j=0; j<ny; j++)
    phi[0][j] = phi[nx][j]; 
}

int main()
{
    int nsteps = 100;
    int nx = 100;
    int ny = 100;
    
    int istep;
    double ** phi = (double **) alloc_2d(nx, ny);
    Parameter params;

    output_tecplot("output.tec", "w", phi, nx, ny, 0);

    for (istep=0; istep<nsteps; istep++)
    {
      solve_step(phi, nx, ny, params);
      pbc(phi, nx, ny);
    }

    free(phi);

    return 0;
}

