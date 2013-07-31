#include <stdlib.h>
#include <stdio.h>

void ** alloc_2d(int nx, int ny)
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

void output_tecplot(double ** phi, int nx, int ny)
{
}

int main()
{
    int nx = 100;
    int ny = 100;

    double ** phi = (double **) alloc_2d(nx, ny);

    free(phi);

    return 0;
}

void array_init ()
{
  void * p;
  int n = 400;
  p = (int *) malloc(n * sizeof(int));
}
