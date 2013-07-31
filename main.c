#include <stdlib.h>
#include <stdio.h>
#include "solve.h"
#include "types.h"

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
}

/* Set initial values of the phi array */
void initial_conditions(double ** phi, int nx, int ny)
{
}

/* Apply Periodic Boundary Conditions */
void pbc(double ** phi, int nx, int ny)
{
}

int main()
{
    int nsteps = 100;
    int nx = 100;
    int ny = 100;

    int istep;
    double ** phi = (double **) alloc_2d(nx, ny);

    for (istep=0; istep<nsteps; istep++)
    {
    }

    free(phi);

    return 0;
}

