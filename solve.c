#include "solve.h"

void solve_step(double ** phi, int nx, int ny, Parameter params)
{
  int i, j;
  double phi1, phi2, phi3, term1, term_xy;
  double dt = params.dt;
  double dx = params.dx;

  for(i=0; i<nx; i++){
    for(j=0; j<ny; j++){
      phi1 = phi[i][j];
      phi2 = phi[i][j] * phi[i][j];
      phi3 = phi[i][j] * phi[i][j] * phi[i][j];

      term1 = 0.1 * (1 - phi2) * (1 - phi2);
      term_xy = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - 4 * phi[i][j]) / (dx * dx);
      phi[i][j] = phi[i][j] + dt * (phi3 - phi1 + term1);
    }
  }
}
