
#include <stdlib.h>
#include <stdio.h>

void ** alloc_2d(int nx, int ny)
{
    int i;
    void ** array = malloc(sizeof(double) * nx);
    for (i=0; i<nx; i++)
        array[i] = malloc(sizeof(double) * ny);

    return array;
}


int main()
{

    int nx = 100;
    int ny = 100;

    double ** phi = (double **) alloc_2d(nx, ny);

    return 0;

}

