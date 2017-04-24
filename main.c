#include <stdio.h>
#include "struct.h"
#include "init.h"
#include "io.h"
#include "diff.h"
#include "converter.h"


int main(int argc, char **argv){
  physics_grid * P_state;
  U_grid * U_state;
  F_grid  * F_state;

  P_state = create_physics_grid();
  U_state = create_U_grid();
  Fx_state = create_F_grid();
  Fy_state = create_F_grid();
  Fz_state = create_F_grid();
  
  init_problem(P_state, U_state, Fx_state,Fy_state,Fz_state, SEDOV);
  initMatrixP(P_state);
  init_UandF(P_state, U_state, Fx_state,Fy_state,Fz_state, SEDOV);
  while()
    {
      VolumenesFinitos(U_state);
    }

  return 0;
}
