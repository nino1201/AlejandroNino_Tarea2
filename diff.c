#include "struct.h"
#include "converter.h"
#include <stdlib.h>
#include <stdio.h>

void CalFx(U_grid *U, F_grid *Fx)
{
  int i;
  int j;
  int k;
  FLOAT N = U->N_cells;

  for(i=0;i<N;i++)
    {
    for(j=0;j<N;j++)
      {
      for(k=0;k<N;k++)
	{
  
	U1=U->U_1[t(i,j,k)];
	U2=U->U_2[t(i,j,k)];
	U3=U->U_3[t(i,j,k)];
	U4=U->U_4[t(i,j,k)];
	U5=U->U_5[t(i,j,k)];
	PR=PR(U1,U2,U3,U4,U5);
  
	Fx->F_1[t(i,j,k)]=U2;
	Fx->F_2[t(i,j,k)]=pow(U2,2)/U1+PR;
	Fx->F_3[t(i,j,k)]=U2*U3/U1;
	Fx->F_4[t(i,j,k)]=U2*U4/U1;
	Fx->F_5[t(i,j,k)]=U2*(U5+PR)/U1;
      }
     }
    }
}

void CalFy( U_grid *U, F_grid *Fy)
{
  int i;
  int j;
  int k;
  FLOAT N = U->N_cells;

  for(i=0;i<N;i++)
    {
    for(j=0;j<N;j++)
      {
      for(k=0;k<N;k++)
	{
  
	U1=U->U_1[t(i,j,k)];
	U2=U->U_2[t(i,j,k)];
	U3=U->U_3[t(i,j,k)];
	U4=U->U_4[t(i,j,k)];
	U5=U->U_5[t(i,j,k)];
	PR=PR(U1,U2,U3,U4,U5);
	      
	Fy->F_1[t(i,j,k)]=U3;
	Fy->F_2[t(i,j,k)]=U2*U3/U1;
	Fy->F_3[t(i,j,k)]=pow(U3,2)/U1+PR;
	Fy->F_4[t(i,j,k)]=U3*U4/U1;
	Fy->F_5[t(i,j,k)]=U3*(U5+PR)/U1;
	
      }
     }
    }
}

void CalFz( U_grid *U, F_grid *Fz)
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
  
	U1=U->U_1[t(i,j,k)];
	U2=U->U_2[t(i,j,k)];
	U3=U->U_3[t(i,j,k)];
	U4=U->U_4[t(i,j,k)];
	U5=U->U_5[t(i,j,k)];
	PR=PR(U1,U2,U3,U4,U5);

	Fz->F_1[t(i,j,k)]=U3;
	Fz->F_2[t(i,j,k)]=U2*U3/U1;
	Fz->F_3[t(i,j,k)]=pow(U3,2)/U1+PR;
	Fz->F_4[t(i,j,k)]=U3*U4/U1;
	Fz->F_5[t(i,j,k)]=U3*(U5+PR)/U1;

	
      }
     }
    }

}
/* Calcula los F */
void CalculateF(physics_grid *P, U_grid *U, F_grid *Fx,F_grid *Fy,F_grid *Fz)
{
  CalFx(U,Fx);
  CalFy(U,Fy);
  CalFz(U,Fz);
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


  CalFx(Ux1,Fx1);
  CalFx(Ux2,Fx2);
  CalFy(Uy1,Fy1);
  CalFy(Uy2,Fy2);
  CalFz(Uz1,Fz1);
  CalFz(Uz2,Fz2);
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
  init_U(Uaux,Paux,SEDOV);
  
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

FLOAT PR(FLOAT U1,FLOAT U2,FLOAT U3,FLOAT U4,FLOAT U5)
{
	return (GAMMA-1)*(U3-0.5*U2*U2/(U1));
}
