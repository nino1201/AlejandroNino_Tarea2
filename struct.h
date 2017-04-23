#ifndef STRUCT_H
#define STRUCT_H

#define GAMMA 1.4

#define SEDOV 1

#define FLOAT double
typedef struct physics_grid_str
{
  FLOAT L_x;
  FLOAT L_y;
  FLOAT L_z;
  FLOAT delta_x;
  FLOAT delta_y;
  FLOAT delta_z;
  int N_x;
  int N_y;
  int N_z;
  int N_cells;
  FLOAT *p;
  FLOAT *rho;
  FLOAT *vx;
  FLOAT *vy;
  FLOAT *vz;
} physics_grid;

typedef struct U_grid_str{
  int N_x;
  int N_y;
  int N_z;
  int N_cells;
  FLOAT *U_1;
  FLOAT *U_2;
  FLOAT *U_3;
  FLOAT *U_4;
  FLOAT *U_5;
} U_grid;


typedef struct G_grid_str{
  int N_x;
  int N_y;
  int N_z;
  int N_cells;
  FLOAT *F_1;
  FLOAT *F_2;
  FLOAT *F_3;
  FLOAT *F_4;
  FLOAT *F_5;

} F_grid;


#endif
