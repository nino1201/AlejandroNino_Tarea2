#include "struct.h"
#include "converter.h"
#include "init.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

FLOAT PR(FLOAT U1,FLOAT U2,FLOAT U3,FLOAT U4,FLOAT U5);
FLOAT time(U_grid *U);
FLOAT cs(FLOAT p,FLOAT rho);
void CalFx(U_grid *U, F_grid *Fx)
{
  int i;
  int N = U->N_cells;

  FLOAT U1;
  FLOAT U2;
  FLOAT U3;
  FLOAT U4;
  FLOAT U5;
  FLOAT pr;
  
  for(i=0;i<N;i++)
  {  
	U1=U->U1[i];
	U2=U->U2[i];
	U3=U->U3[i];
	U4=U->U4[i];
	U5=U->U5[i];
	pr=PR(U1,U2,U3,U4,U5);
  
	Fx->F1[i]=U2;
	Fx->F2[i]=pow(U2,2)/U1+pr;
	Fx->F3[i]=U2*U3/U1;
	Fx->F4[i]=U2*U4/U1;
	Fx->F5[i]=U2*(U5+pr)/U1;
  }
}

void CalFy( U_grid *U, F_grid *Fy)
{
  int i;
  int N = U->N_cells;

  FLOAT U1;
  FLOAT U2;
  FLOAT U3;
  FLOAT U4;
  FLOAT U5;
  FLOAT pr;
  
  
  for(i=0;i<N;i++)
  {  
	U1=U->U1[i];
	U2=U->U2[i];
	U3=U->U3[i];
	U4=U->U4[i];
	U5=U->U5[i];
	pr=PR(U1,U2,U3,U4,U5);
	      
	Fy->F1[i]=U3;
	Fy->F2[i]=U2*U3/U1;
	Fy->F3[i]=pow(U3,2)/U1+pr;
	Fy->F4[i]=U3*U4/U1;
	Fy->F5[i]=U3*(U5+pr)/U1;
  }
}

void CalFz( U_grid *U, F_grid *Fz)
{
  int i;
  FLOAT U1;
  FLOAT U2;
  FLOAT U3;
  FLOAT U4;
  FLOAT U5;
  FLOAT pr;
  
  int N = U->N_cells;

  for(i=0;i<N;i++)
  {
	U1=U->U1[i];
	U2=U->U2[i];
	U3=U->U3[i];
	U4=U->U4[i];
	U5=U->U5[i];
	pr=PR(U1,U2,U3,U4,U5);

	Fz->F1[i]=U3;
	Fz->F2[i]=U2*U3/U1;
	Fz->F3[i]=pow(U3,2)/U1+pr;
	Fz->F4[i]=U3*U4/U1;
	Fz->F5[i]=U3*(U5+pr)/U1;

    }
}
/* Calcula los F */
void CalculateF(U_grid *U, F_grid *Fx,F_grid *Fy,F_grid *Fz)
{
  CalFx(U,Fx);
  CalFy(U,Fy);
  CalFz(U,Fz);
}

/* Calcula todos los F gorritos */
void CalculateFG(physics_grid *P, U_grid *U,F_grid *Fx1,F_grid *Fx2,F_grid *Fy1,F_grid *Fy2,F_grid *Fz1,F_grid *Fz2 )
{
  U_grid *Ux1;
  U_grid *Ux2;
  U_grid *Uy1;
  U_grid *Uy2;
  U_grid *Uz1;
  U_grid *Uz2;
  
  Ux1 = create_U_grid();
  Ux2 = create_U_grid();
  Uy1 = create_U_grid();
  Uy2 = create_U_grid();
  Uz1 = create_U_grid();
  Uz2 = create_U_grid();
  
  init_U(Ux1,P,SEDOV);
  init_U(Ux2,P,SEDOV);
  init_U(Uy1,P,SEDOV);
  init_U(Uy2,P,SEDOV);
  init_U(Uz1,P,SEDOV);
  init_U(Uz2,P,SEDOV);

  init_Uc(P,Ux1,SEDOV);
  init_Uc(P,Ux2,SEDOV);
  init_Uc(P,Uy1,SEDOV);
  init_Uc(P,Uy2,SEDOV);
  init_Uc(P,Uz1,SEDOV);
  init_Uc(P,Uz2,SEDOV);
  
  
  int i;
  int j;
  int k;

  int N=U->N_x;
  for(i=1;i<N-1;i++){
    for(j=1;j<N-1;j++){
      for(k=1;k<N-1;k++){
	Ux1->U1[t(i,j,k)]=(U->U1[t(i-1,j,k)]+U->U1[t(i,j,k)])/2;
	Uy1->U1[t(i,j,k)]=(U->U1[t(i,j-1,k)]+U->U1[t(i,j,k)])/2;
	Uz1->U1[t(i,j,k)]=(U->U1[t(i,j,k-1)]+U->U1[t(i,j,k)])/2;
	Ux2->U1[t(i,j,k)]=(U->U1[t(i+1,j,k)]+U->U1[t(i,j,k)])/2;
	Uy2->U1[t(i,j,k)]=(U->U1[t(i,j+1,k)]+U->U1[t(i,j,k)])/2;
	Uz2->U1[t(i,j,k)]=(U->U1[t(i,j,k+1)]+U->U1[t(i,j,k)])/2;

	Ux1->U2[t(i,j,k)]=(U->U2[t(i-1,j,k)]+U->U2[t(i,j,k)])/2;
	Uy1->U2[t(i,j,k)]=(U->U2[t(i,j-1,k)]+U->U2[t(i,j,k)])/2;
	Uz1->U2[t(i,j,k)]=(U->U2[t(i,j,k-1)]+U->U2[t(i,j,k)])/2;
	Ux2->U2[t(i,j,k)]=(U->U2[t(i+1,j,k)]+U->U2[t(i,j,k)])/2;
	Uy2->U2[t(i,j,k)]=(U->U2[t(i,j+1,k)]+U->U2[t(i,j,k)])/2;
	Uz2->U2[t(i,j,k)]=(U->U2[t(i,j,k+1)]+U->U2[t(i,j,k)])/2;

	Ux1->U3[t(i,j,k)]=(U->U3[t(i-1,j,k)]+U->U3[t(i,j,k)])/2;
	Uy1->U3[t(i,j,k)]=(U->U3[t(i,j-1,k)]+U->U3[t(i,j,k)])/2;
	Uz1->U3[t(i,j,k)]=(U->U3[t(i,j,k-1)]+U->U3[t(i,j,k)])/2;
	Ux2->U3[t(i,j,k)]=(U->U3[t(i+1,j,k)]+U->U3[t(i,j,k)])/2;
	Uy2->U3[t(i,j,k)]=(U->U3[t(i,j+1,k)]+U->U3[t(i,j,k)])/2;
	Uz2->U3[t(i,j,k)]=(U->U3[t(i,j,k+1)]+U->U3[t(i,j,k)])/2;

	Ux1->U4[t(i,j,k)]=(U->U4[t(i-1,j,k)]+U->U4[t(i,j,k)])/2;
	Uy1->U4[t(i,j,k)]=(U->U4[t(i,j-1,k)]+U->U4[t(i,j,k)])/2;
	Uz1->U4[t(i,j,k)]=(U->U4[t(i,j,k-1)]+U->U4[t(i,j,k)])/2;
	Ux2->U4[t(i,j,k)]=(U->U4[t(i+1,j,k)]+U->U4[t(i,j,k)])/2;
	Uy2->U4[t(i,j,k)]=(U->U4[t(i,j+1,k)]+U->U4[t(i,j,k)])/2;
	Uz2->U4[t(i,j,k)]=(U->U4[t(i,j,k+1)]+U->U4[t(i,j,k)])/2;

	Ux1->U5[t(i,j,k)]=(U->U5[t(i-1,j,k)]+U->U5[t(i,j,k)])/2;
	Uy1->U5[t(i,j,k)]=(U->U5[t(i,j-1,k)]+U->U5[t(i,j,k)])/2;
	Uz1->U5[t(i,j,k)]=(U->U5[t(i,j,k-1)]+U->U5[t(i,j,k)])/2;
	Ux2->U5[t(i,j,k)]=(U->U5[t(i+1,j,k)]+U->U5[t(i,j,k)])/2;
	Uy2->U5[t(i,j,k)]=(U->U5[t(i,j+1,k)]+U->U5[t(i,j,k)])/2;
	Uz2->U5[t(i,j,k)]=(U->U5[t(i,j,k+1)]+U->U5[t(i,j,k)])/2;
	
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
  physics_grid * Paux;
  U_grid * U_aux;
  F_grid  * Fx1;
  F_grid  * Fx2;
  F_grid  * Fy1;
  F_grid  * Fy2;
  F_grid  * Fz1;
  F_grid  * Fz2;

  Paux = create_physics_grid();
  init_P(Paux,SEDOV);
  U_aux = create_U_grid();
  init_U(U_aux,Paux,SEDOV);
  
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

  /* Calcula todos los F gorrito */
  CalculateFG(Paux,U,Fx1,Fx2,Fy1,Fy2,Fz1,Fz2);
  
  int i;
  int Nc=U->N_cells;
  FLOAT dt=time(U);  
  for(i=0;i<Nc;i++)
    {
      U_aux->U1[i]= U->U1[i];
      U_aux->U2[i]= U->U2[i];
      U_aux->U3[i]= U->U3[i];
      U_aux->U4[i]= U->U4[i];
      U_aux->U5[i]= U->U5[i];
    }

  for(i=0;i<Nc;i++)
  {
	U->U1[i]=U_aux->U1[i]+dt/dx*(Fx1->F1[i]-Fx2->F1[i])+dt/dx*(Fy1->F1[i]-Fy2->F1[i])+dt/dx*(Fz1->F1[i]-Fz2->F1[i]);
	U->U2[i]=U_aux->U2[i]+dt/dx*(Fx1->F2[i]-Fx2->F2[i])+dt/dx*(Fy1->F2[i]-Fy2->F2[i])+dt/dx*(Fz1->F2[i]-Fz2->F2[i]);
	U->U3[i]=U_aux->U3[i]+dt/dx*(Fx1->F3[i]-Fx2->F3[i])+dt/dx*(Fy1->F3[i]-Fy2->F3[i])+dt/dx*(Fz1->F3[i]-Fz2->F3[i]);
	U->U4[i]=U_aux->U4[i]+dt/dx*(Fx1->F4[i]-Fx2->F4[i])+dt/dx*(Fy1->F4[i]-Fy2->F4[i])+dt/dx*(Fz1->F4[i]-Fz2->F4[i]);
	U->U5[i]=U_aux->U5[i]+dt/dx*(Fx1->F5[i]-Fx2->F5[i])+dt/dx*(Fy1->F5[i]-Fy2->F5[i])+dt/dx*(Fz1->F5[i]-Fz2->F5[i]);

	/* halla la velocidad maxima y la guarda */
	if(U->vx_MAX<(U->U2[i]/U->U1[i]))
	  {
	    U->vx_MAX=U->U2[i]/U->U1[i];
	    U->rhox_MAX=U->U1[i];
	    U->px_MAX=PR(U->U1[i],U->U2[i],U->U3[i],U->U4[i],U->U5[i]);
	  }
	if(U->vy_MAX<(U->U2[i]/U->U1[i]))
	  {
	    U->vy_MAX=U->U2[i]/U->U1[i];
	    U->rhox_MAX=U->U1[i];
	    U->px_MAX=PR(U->U1[i],U->U2[i],U->U3[i],U->U4[i],U->U5[i]);
	  }
	if(U->vy_MAX<(U->U2[i]/U->U1[i]))
	  {
	    U->vy_MAX=U->U2[i]/U->U1[i];
	    U->rhox_MAX=U->U1[i];
	    U->px_MAX=PR(U->U1[i],U->U2[i],U->U3[i],U->U4[i],U->U5[i]);
	  }
  }
  /* intento velocidades maximas */
  
}

/* velocidad de la luz */
FLOAT cs(FLOAT p,FLOAT rho)
{
  return pow(GAMMA*(p/rho),0.5);
}

/* debe retornar el dt para velocidades maximas */
FLOAT time(U_grid *U) 
 { 
   FLOAT max;
   max=U->vx_MAX+cs(U->px_MAX,U->rhox_MAX);
   if(max<U->vy_MAX+cs(U->py_MAX,U->rhoy_MAX))
     {
       max=U->vy_MAX+cs(U->py_MAX,U->rhoy_MAX);
     }
   if(max<U->vz_MAX+cs(U->pz_MAX,U->rhoz_MAX))
     {
       max=U->vz_MAX+cs(U->pz_MAX,U->rhoz_MAX);
     }
   return dx/max; 
 } 

 /* calcula la energia */

FLOAT E(FLOAT rho,FLOAT e,FLOAT u,FLOAT v,FLOAT w)
{
  return rho*e+(pow(u,2)+pow(v,2)+pow(w,2))/2;
}

/* calcula la enrgia interna */

FLOAT e(FLOAT p,FLOAT rho)
{
  return p/(rho*(GAMMA-1));
}

FLOAT PR(FLOAT U1,FLOAT U2,FLOAT U3,FLOAT U4,FLOAT U5)
{
	return (GAMMA-1)*(U3-0.5*U2*U2/(U1));
}
