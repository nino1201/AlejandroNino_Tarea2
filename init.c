#include "struct.h"
#include "converter.h"
#include <stdlib.h>
#include <stdio.h>

void init_to_zero(FLOAT *p, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    p[i] = 0.0;
  }
}


physics_grid * create_physics_grid(void){
  physics_grid *G;
  if(!(G = malloc(sizeof(physics_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  } 
  
  G->L_x=0.0;
  G->L_y=0.0;
  G->L_z=0.0;
  G->delta_x=0.0;
  G->delta_y=0.0;
  G->delta_z=0.0;
  G->N_x=0;
  G->N_y=0;
  G->N_z=0;
  G->N_cells=0;
  G->p=NULL;
  G->rho=NULL;
  G->vx=NULL;
  G->vy=NULL;
  G->vz=NULL;
  return G;
}


U_grid * create_U_grid(void){
  U_grid *G;
  if(!(G = malloc(sizeof(U_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  } 
  G->N_x=0.0;
  G->N_y=0.0;
  G->N_z=0.0;
  G->N_cells=0.0;
  G->U_1=NULL;
  G->U_2=NULL;
  G->U_3=NULL;
  G->U_4=NULL;
  G->U_5=NULL;
  return G;
}

F_grid * create_F_grid(void){
  F_grid *G;
  if(!(G = malloc(sizeof(F_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  } 
  G->N_x=0.0;
  G->N_y=0.0;
  G->N_z=0.0;
  G->N_cells=0.0;
  G->F_1=NULL;
  G->F_2=NULL;
  G->F_3=NULL;
  G->F_4=NULL;
  G->F_5=NULL;

  return G;
}


void init_P(physics_grid *P, int problem){
  
  P->L_x = 250.0;
  P->L_y = 250.0;
  P->L_z = 250.0;    
  P->delta_x = 5.0;
  P->delta_y = 5.0;
  P->delta_z = 5.0;
  P->N_x = (int)(P->L_x/P->delta_x);
  P->N_y = (int)(P->L_y/P->delta_y);
  P->N_z = (int)(P->L_z/P->delta_z);
  P->N_cells = P->N_x * P->N_y * P->N_z;

  if(!(P->p=malloc(P->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  init_to_zero(P->p, P->N_cells);
  
  if(!(P->rho=malloc(P->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  init_to_zero(P->rho, P->N_cells);
  
  if(!(P->vx=malloc(P->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with vx allocation");
    exit(1);
  }
  init_to_zero(P->vx, P->N_cells);
  
  if(!(P->vy=malloc(P->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with vy allocation");
    exit(1);
  }
  init_to_zero(P->vy, P->N_cells);
  
  if(!(P->vz=malloc(P->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with vz allocation");
    exit(1);
  }    
  init_to_zero(P->vz, P->N_cells);
}

void init_U(U_grid *U,physics_grid *P,int problem)
{
  
  U->N_x = P->N_x;
  U->N_y = P->N_y;
  U->N_z = P->N_z;
  U->N_cells = P->N_cells;

 
  if(!(U->U_1=malloc(U->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_1, U->N_cells);
  
  if(!(U->U_2=malloc(U->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_2 allocation");
    exit(1);
  }
  init_to_zero(U->U_2, U->N_cells);
  
  if(!(U->U_3=malloc(U->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_2 allocation");
    exit(1);
  }
  init_to_zero(U->U_3, U->N_cells);
  
  if(!(U->U_4=malloc(U->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_2 allocation");
    exit(1);
  }
  init_to_zero(U->U_4, U->N_cells);
  
  if(!(U->U_5=malloc(U->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_2 allocation");
    exit(1);
  }
  init_to_zero(U->U_5, U->N_cells);

}

void init_F(F_grid *F,physics_grid *P,int problem){
  F->N_x = P->N_x;
  F->N_y = P->N_y;
  F->N_z = P->N_z;
  F->N_cells = P->N_cells;
  
  
  if(!(F->F_1=malloc(F->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_1_X allocation");
    exit(1);
  }
  init_to_zero(F->F_1, F->N_cells);
  
  if(!(F->F_2=malloc(F->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_2_X allocation");
    exit(1);
  }
  init_to_zero(F->F_2, F->N_cells);
  
  if(!(F->F_3=malloc(F->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_3_X allocation");
    exit(1);
  }
  init_to_zero(F->F_3, F->N_cells);
  
  if(!(F->F_4=malloc(F->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_4_X allocation");
    exit(1);
  }
  init_to_zero(F->F_4, F->N_cells);
  
  if(!(F->F_5=malloc(F->N_cells * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_5_X allocation");
    exit(1);
  }
  init_to_zero(F->F_5, F->N_cells);
  
}
