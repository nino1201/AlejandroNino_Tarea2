#include <stdio.h>
#include "struct.h"
#include "init.h"
#include "io.h"
#include "diff.h"
#include "converter.h"


int main(int argc, char **argv){
  physics_grid * P_state;
  U_grid * U_state;
  F_grid  * Fx_state;
  F_grid  * Fy_state;
  F_grid  * Fz_state;
  
  P_state = create_physics_grid();
  U_state = create_U_grid();
  Fx_state = create_F_grid();
  Fy_state = create_F_grid();
  Fz_state = create_F_grid();
  
  init_problem(P_state, U_state, Fx_state,Fy_state,Fz_state, SEDOV);
  initMatrixP(P_state,SEDOV);
  init_Uc(P_state, U_state,SEDOV);
  init_Fc(P_state,Fx_state,Fy_state,Fz_state,SEDOV);
  int j;
  for(j=0;j<2;j++)
  {
    VolumenesFinitos(U_state);
  }
  for(j=0;j<U_state->N_cells;j++)
  {
    printf("%f %f %f %f %f\n",U_state->U1[j],U_state->U2[j],U_state->U3[j],U_state->U4[j],U_state->U5[j]); 
  }
  return 0;
}
