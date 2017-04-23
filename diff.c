#include "struct.h"
#include "converter.h"
#include <stdlib.h>
#include <stdio.h>

void CalFx(physics_grid *P, U_grid *U, F_grid *Fx)
{
  int i;
  int j;
  int k;
  FLOAT N = P->N_cells;

  for(i=0;i<N;i++)
    {
    for(j=0;j<N;j++)
      {
      for(k=0;k<N;k++)
	{
  
	FLOAT U1=U->U_1[t(i,j,k)];
	FLOAT U2=U->U_2[t(i,j,k)];
	FLOAT U3=U->U_3[t(i,j,k)];
	FLOAT U4=U->U_4[t(i,j,k)];
	FLOAT U5=U->U_5[t(i,j,k)];
	FLOAT PR=P->p[t(i,j,k)];
  
	FLOAT Fx->F_1[t(i,j,k)]=U2;
	FLOAT Fx->F_2[t(i,j,k)]=pow(U2,2)/U1+PR;
	FLOAT Fx->F_3[t(i,j,k)]=U2*U3/U1;
	FLOAT Fx->F_4[t(i,j,k)]=U2*U4/U1;
	FLOAT Fx->F_5[t(i,j,k)]=U2*(U5+PR)/U1;
      }
     }
    }
}

void CalFy(physics_grid *P, U_grid *U, F_grid *Fy)
{
  int i;
  int j;
  int k;
  FLOAT N = P->N_cells;

  for(i=0;i<N;i++)
    {
    for(j=0;j<N;j++)
      {
      for(k=0;k<N;k++)
	{
  
	FLOAT U1=U->U_1[t(i,j,k)];
	FLOAT U2=U->U_2[t(i,j,k)];
	FLOAT U3=U->U_3[t(i,j,k)];
	FLOAT U4=U->U_4[t(i,j,k)];
	FLOAT U5=U->U_5[t(i,j,k)];
	FLOAT PR=P->p[t(i,j,k)];

	FLOAT Fy->F_1[t(i,j,k)]=U3;
	FLOAT Fy->F_2[t(i,j,k)]=U2*U3/U1;
	FLOAT Fy->F_3[t(i,j,k)]=pow(U3,2)/U1+PR;
	FLOAT Fy->F_4[t(i,j,k)]=U3*U4/U1;
	FLOAT Fy->F_5[t(i,j,k)]=U3*(U5+PR)/U1;

	
      }
     }
    }
}

void CalFz(physics_grid *P, U_grid *U, F_grid *Fz)
{
  int i;
  int j;
  int k;
  FLOAT N = P->N_cells;

  for(i=0;i<N;i++)
    {
    for(j=0;j<N;j++)
      {
      for(k=0;k<N;k++)
	{
  
	FLOAT U1=U->U_1[t(i,j,k)];
	FLOAT U2=U->U_2[t(i,j,k)];
	FLOAT U3=U->U_3[t(i,j,k)];
	FLOAT U4=U->U_4[t(i,j,k)];
	FLOAT U5=U->U_5[t(i,j,k)];
	FLOAT PR=P->p[t(i,j,k)];

	FLOAT Fz->F_1[t(i,j,k)]=U3;
	FLOAT Fz->F_2[t(i,j,k)]=U2*U3/U1;
	FLOAT Fz->F_3[t(i,j,k)]=pow(U3,2)/U1+PR;
	FLOAT Fz->F_4[t(i,j,k)]=U3*U4/U1;
	FLOAT Fz->F_5[t(i,j,k)]=U3*(U5+PR)/U1;

	
      }
     }
    }

}
/* Calcula los F */
void CalculateF(physics_grid *P, U_grid *U, F_grid *Fx,F_grid *Fy,F_grid *Fz)
{
  CalFx(P,U,Fx);
  CalFy(P,U,Fy);
  CalFz(P,U,Fz);
}

/* Calcula todos los F gorritos */
void CalculateFG(physics_grid *P, U_grid *U,F_grid *Fx1,F_grid *Fx2,F_grid *Fy1,F_grid *Fy2,F_grid *Fz1,F_grid *Fz2 )
{
  U_gird *Ux1;
  U_gird *Ux2;
  U_gird *Uy1;
  U_gird *Uy2;
  U_gird *Uz1;
  U_gird *Uz2;
  
  Ux1 = create_U_grid();
  Ux2 = create_U_grid();
  Uy1 = create_U_grid();
  Uy2 = create_U_grid();
  Uz1 = create_U_grid();
  Uz2 = create_U_grid();
  
  init_U(P,Ux1,SEDOV);
  init_U(P,Ux2,SEDOV);
  init_U(P,Uy1,SEDOV);
  init_U(P,Uy2,SEDOV);
  init_U(P,Uz1,SEDOV);
  init_U(P,Uz2,SEDOV);
  
  
  int i;
  int j;
  int k;

  
  for(i=1;i<N-1;i++){
    for(j=1;j<N-1;j++){
      for(k=1;k<N-1;k++){
	Ux1->U_1[t(i,j,k)]=(U->U_1[i-1,j,k]+U->U_1[i,j,k])/2;
	Uy1->U_1[t(i,j,k)]=(U->U_1[i,j-1,k]+U->U_1[i,j,k])/2;
	Uz1->U_1[t(i,j,k)]=(U->U_1[i,j,k-1]+U->U_1[i,j,k])/2;
	Ux2->U_1[t(i,j,k)]=(U->U_1[i+1,j,k]+U->U_1[i,j,k])/2;
	Uy2->U_1[t(i,j,k)]=(U->U_1[i,j+1,k]+U->U_1[i,j,k])/2;
	Uz2->U_1[t(i,j,k)]=(U->U_1[i,j,k+1]+U->U_1[i,j,k])/2;

	Ux1->U_2[t(i,j,k)]=(U->U_2[i-1,j,k]+U->U_2[i,j,k])/2;
	Uy1->U_2[t(i,j,k)]=(U->U_2[i,j-1,k]+U->U_2[i,j,k])/2;
	Uz1->U_2[t(i,j,k)]=(U->U_2[i,j,k-1]+U->U_2[i,j,k])/2;
	Ux2->U_2[t(i,j,k)]=(U->U_2[i+1,j,k]+U->U_2[i,j,k])/2;
	Uy2->U_2[t(i,j,k)]=(U->U_2[i,j+1,k]+U->U_2[i,j,k])/2;
	Uz2->U_2[t(i,j,k)]=(U->U_2[i,j,k+1]+U->U_2[i,j,k])/2;

	Ux1->U_3[t(i,j,k)]=(U->U_3[i-1,j,k]+U->U_3[i,j,k])/2;
	Uy1->U_3[t(i,j,k)]=(U->U_3[i,j-1,k]+U->U_3[i,j,k])/2;
	Uz1->U_3[t(i,j,k)]=(U->U_3[i,j,k-1]+U->U_3[i,j,k])/2;
	Ux2->U_3[t(i,j,k)]=(U->U_3[i+1,j,k]+U->U_3[i,j,k])/2;
	Uy2->U_3[t(i,j,k)]=(U->U_3[i,j+1,k]+U->U_3[i,j,k])/2;
	Uz2->U_3[t(i,j,k)]=(U->U_3[i,j,k+1]+U->U_3[i,j,k])/2;

	Ux1->U_4[t(i,j,k)]=(U->U_4[i-1,j,k]+U->U_4[i,j,k])/2;
	Uy1->U_4[t(i,j,k)]=(U->U_4[i,j-1,k]+U->U_4[i,j,k])/2;
	Uz1->U_4[t(i,j,k)]=(U->U_4[i,j,k-1]+U->U_4[i,j,k])/2;
	Ux2->U_4[t(i,j,k)]=(U->U_4[i+1,j,k]+U->U_4[i,j,k])/2;
	Uy2->U_4[t(i,j,k)]=(U->U_4[i,j+1,k]+U->U_4[i,j,k])/2;
	Uz2->U_4[t(i,j,k)]=(U->U_4[i,j,k+1]+U->U_4[i,j,k])/2;

	Ux1->U_5[t(i,j,k)]=(U->U_5[i-1,j,k]+U->U_5[i,j,k])/2;
	Uy1->U_5[t(i,j,k)]=(U->U_5[i,j-1,k]+U->U_5[i,j,k])/2;
	Uz1->U_5[t(i,j,k)]=(U->U_5[i,j,k-1]+U->U_5[i,j,k])/2;
	Ux2->U_5[t(i,j,k)]=(U->U_5[i+1,j,k]+U->U_5[i,j,k])/2;
	Uy2->U_5[t(i,j,k)]=(U->U_5[i,j+1,k]+U->U_5[i,j,k])/2;
	Uz2->U_5[t(i,j,k)]=(U->U_5[i,j,k+1]+U->U_5[i,j,k])/2;
	
      }
    }
  }

  CalFx(P,Ux1,Fx1);
  CalFx(P,Ux2,Fx2);
  CalFy(P,Uy1,Fy1);
  CalFy(P,Uy2,Fy2);
  CalFz(P,Uz1,Fz1);
  CalFz(P,Uz2,Fz2);
}

void VolumenesFinitos( U_grid *U)
{
  physics_grid * P_aux;
  U_grid * U_aux;
  F_grid  * Fx1;
  F_grid  * Fx2;
  F_grid  * Fy1;
  F_grid  * Fy2;
  F_grid  * Fz1;
  F_grid  * Fz2;

  P_aux = create_physics_grid();
  init_P(Paux,SEDOV);
  U_aux = create_U_grid();
  init_U(Uaux,SEDOV);
  
  Fx1 = create_F_grid();
  Fx2 = create_F_grid();
  Fy1 = create_F_grid();
  Fy2 = create_F_grid();
  Fz1 = create_F_grid();
  Fz2 = create_F_grid();

  init_F(Fx1,Paux,SEDOV);
  init_F(Fx2,Paux,SEDOV);
  init_F(Fy1,Paux,SEDOV);
  init_F(Fy2,Paux,SEDOV);
  init_F(Fz1,Paux,SEDOV);
  init_F(Fz2,Paux,SEDOV);
  
  CalculateFG(U,Fx1,Fx2,Fy1,Fy2,Fz1,Fz2);
  
  int i;
  int j;
  int k;

  FLOAT dt=time(U);  
  for(i=0;i<N;i++)
    {
      for(i=0;i<N;i++)
	{
	  for(i=0;i<N;i++)
	    {
	      U_aux->U_1[t(i,j,k)]= U->U_1[t(i,j,k)];
	      U_aux->U_2[t(i,j,k)]= U->U_2[t(i,j,k)];
	      U_aux->U_3[t(i,j,k)]= U->U_3[t(i,j,k)];
	      U_aux->U_4[t(i,j,k)]= U->U_4[t(i,j,k)];
	      U_aux->U_5[t(i,j,k)]= U->U_5[t(i,j,k)];
	    }
	} 
    }

  
  for(i=0;i<N;i++)
    {
      for(i=0;i<N;i++)
	{
	  for(i=0;i<N;i++)
	    {

	      U->U_1[t(i,j,k)]=U_aux->U_1[t(i,j,k)]+dt/dx*(Fx1->F_1[t(i,j,k)]-Fx2->F_1[t(i,j,k)])+dt/dx*(Fy1->F_1[t(i,j,k)]-Fy2->F_1[t(i,j,k)])+dt/dx*(Fz1->F_1[t(i,j,k)]-Fz2->F_1[t(i,j,k)]);
	      U->U_2[t(i,j,k)]=U_aux->U_2[t(i,j,k)]+dt/dx*(Fx1->F_2[t(i,j,k)]-Fx2->F_2[t(i,j,k)])+dt/dx*(Fy1->F_2[t(i,j,k)]-Fy2->F_2[t(i,j,k)])+dt/dx*(Fz1->F_2[t(i,j,k)]-Fz2->F_2[t(i,j,k)]);
	      U->U_3[t(i,j,k)]=U_aux->U_3[t(i,j,k)]+dt/dx*(Fx1->F_3[t(i,j,k)]-Fx2->F_3[t(i,j,k)])+dt/dx*(Fy1->F_3[t(i,j,k)]-Fy2->F_3[t(i,j,k)])+dt/dx*(Fz1->F_3[t(i,j,k)]-Fz2->F_3[t(i,j,k)]);
	      U->U_4[t(i,j,k)]=U_aux->U_4[t(i,j,k)]+dt/dx*(Fx1->F_4[t(i,j,k)]-Fx2->F_4[t(i,j,k)])+dt/dx*(Fy1->F_4[t(i,j,k)]-Fy2->F_4[t(i,j,k)])+dt/dx*(Fz1->F_4[t(i,j,k)]-Fz2->F_4[t(i,j,k)]);
	      U->U_5[t(i,j,k)]=U_aux->U_5[t(i,j,k)]+dt/dx*(Fx1->F_5[t(i,j,k)]-Fx2->F_5[t(i,j,k)])+dt/dx*(Fy1->F_5[t(i,j,k)]-Fy2->F_5[t(i,j,k)])+dt/dx*(Fz1->F_5[t(i,j,k)]-Fz2->F_5[t(i,j,k)])
	      
	    }
	}
    }
}

/* debe retornar el dt para velocidades maximas */
FLOAT time(U_grid *U)
{
  return 0;
}

/* calcula la energia */

FLOAT E(FLOAT rho,FLOAT e,FLOAT u,FLOAT v,FLOAT w)
{
  return rho*e+(pow(u,2)+pow(v,2)+pow(w,2))/2
}

/* calcula la enrgia interna */

FLOAT e(FLOAT p,FLOAT rho)
{
  return p/(rho*(GAMMA-1))
}

