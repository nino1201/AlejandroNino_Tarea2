#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TYPE double
#define dx 0.01


void inicializar(TYPE *p, TYPE *rho, TYPE *v, TYPE *x, int N);
void inicialUF(TYPE *u1,TYPE *u2,TYPE *u3,TYPE *f1,TYPE *f2,TYPE *f3, TYPE *rho, TYPE *v, TYPE *p, int N);
TYPE e(TYPE rho,TYPE p);
TYPE E(TYPE v,TYPE rho,TYPE e);
void LaxWendroff(TYPE *u1,TYPE *u2,TYPE *u3,TYPE *f1,TYPE *f2,TYPE *f3, int N);
TYPE P(TYPE u1,TYPE u2,TYPE u3);
void actualizaF(TYPE *u1,TYPE *u2,TYPE *u3,TYPE *f1,TYPE *f2,TYPE *f3, int N);



int main()
{
  int l = 1;
  int N = l/dx;
  
  TYPE *p= malloc(N*sizeof(TYPE));
  TYPE *rho= malloc(N*sizeof(TYPE));
  TYPE *v= malloc(N*sizeof(TYPE));
  TYPE *x= malloc(N*sizeof(TYPE));

  inicializar(p,rho,v,x,N);

  TYPE *u1= malloc(N*sizeof(TYPE));
  TYPE *u2= malloc(N*sizeof(TYPE));
  TYPE *u3= malloc(N*sizeof(TYPE));
  TYPE *f1= malloc(N*sizeof(TYPE));
  TYPE *f2= malloc(N*sizeof(TYPE));
  TYPE *f3= malloc(N*sizeof(TYPE));

  inicialUF(u1,u2,u3,f1,f2,f3,rho,v,p,N);
  int j=0;
  /* Falta poner la condicion para de que llegue a x=0.9 */
  while(/*abs((u2[(int)(0.9/dx)]/u1[(int)(0.9/dx)])-(u2[(int)(0.9/dx+1)]/u1[(int)(0.9/dx+1)]))>1E-5 &&*/ p[(int)(0.9/dx)]!=1)
    /* en la condicion del while se compara la velocidad en el punto x=0.9 con el siguiente punto*/
    /*basicamente las velocidades son constantes y si existe una discontinuidad va a parar el ciclo */
    {
      LaxWendroff(u1,u2,u3,f1,f2,f3,N);
      actualizaF(u1,u2,u3,f1,f2,f3,N);
      j=j+1;
      printf("%f\n",j);
    }
  int i;
  for(i=0;i<N;i++)
   {
     /* Imprime posicion, densidad, presion y velocidad */
     printf("%f %f %f %f\n",x[i],u1[i],P(u1[i],u2[i],u3[i]),u2[i]/u1[i]);
   }
  return 0;
  
}
/* Inicializa densidad, presion, velocidad y la  posicion  */

void inicializar(TYPE *p,TYPE *rho, TYPE *v,TYPE *x, int N)
{
  int Naux=0.5/dx;
  int i;
  for(i=0;i<N;i++)
    {
      if(i<=Naux)
	{
	  p[i]=1.0;
	  rho[i]=1.0;
	}
      else
	{
	  p[i]=0.1;
	  rho[i]=0.1;
	}
      x[i]=i*dx;
      v[i]=0;
    }
}

/* Inicializa U y F */

void inicialUF(TYPE *u1,TYPE *u2,TYPE *u3,TYPE *f1,TYPE *f2,TYPE *f3, TYPE *rho, TYPE *v, TYPE *p, int N)
{
  int i;
  for(i=0;i<N;i++)
    {
      u1[i]=rho[i];
      u2[i]=rho[i]*v[i];
      u3[i]=E(v[i],rho[i],e(rho[i],p[i]));

      f1[i]=rho[i]*v[i];
      f2[i]=p[i]+rho[i]*pow(v[i],2);
      f3[i]=v[i]*(E(v[i],rho[i],e(rho[i],p[i]))+p[i]);
    } 
}

/* Calcula energia interna */

TYPE e(TYPE rho,TYPE p)
{
  TYPE g=0.4;
  return p/(rho*g);
}

/*Calcula Energia */

TYPE E(TYPE v,TYPE rho,TYPE e)
{
  return rho*e+rho*(pow(v,2))/2;
}

/* realiza el metodo lax-wedroff */

void LaxWendroff(TYPE *u1,TYPE *u2,TYPE *u3,TYPE *f1,TYPE *f2,TYPE *f3, int N)
{
  TYPE* tempu1= malloc(N*sizeof(TYPE));
  TYPE* tempu2= malloc(N*sizeof(TYPE));
  TYPE* tempu3= malloc(N*sizeof(TYPE));

  tempu1[0]=u1[0];
  tempu2[0]=u2[0];
  tempu3[0]=u3[0];

  tempu1[(int)(N-1)]=u1[(int)(N-1)];
  tempu2[(int)(N-1)]=u2[(int)(N-1)];
  tempu3[(int)(N-1)]=u3[(int)(N-1)];
 
  TYPE umax=0;
  TYPE dt=0;
  int i;
  for (i=1;i<N-1;i++)
    {
      if(u1[i]>=u2[i] && u1[i]>=u3[i])
	{
	  umax=u1[i];
	}
      else if(u2[i]>=u1[i] && u2[i]>=u3[i])
	{
	  umax=u2[i];
	}
      else
	{
	  umax=u3[i];
	}
      dt=dx/umax;
      if(dt==1)
	{
	  dt=dt-0.1;
	}
      tempu1[i]=u1[i]-(0.5*dt/dx)*(f1[i+1]-f1[i-1])+(0.25*pow(dt/dx,2))*((u1[i+1]+u1[i])*(f1[i+1]-f1[i])-(u1[i]+u1[i-1])*(f1[i]-f1[i-1]));
      tempu2[i]=u2[i]-(0.5*dt/dx)*(f2[i+1]-f2[i-1])+(0.25*pow(dt/dx,2))*((u2[i+1]+u2[i])*(f2[i+1]-f2[i])-(u2[i]+u2[i-1])*(f2[i]-f2[i-1]));
      tempu3[i]=u3[i]-(0.5*dt/dx)*(f3[i+1]-f3[i-1])+(0.25*pow(dt/dx,2))*((u3[i+1]+u3[i])*(f3[i+1]-f3[i])-(u3[i]+u3[i-1])*(f3[i]-f3[i-1]));
    }
  u1=tempu1;
  u2=tempu2;
  u3=tempu3;

  /* Falta mirar la condicion dt/dx umax<1 */
}


/* Halla la presion */

TYPE P(TYPE u1,TYPE u2,TYPE u3)
{
  return 0.4*(u3-pow(u2,2)/(2.*u1));
}

/* Actualiza el vector F */

void actualizaF(TYPE *u1,TYPE *u2,TYPE *u3,TYPE *f1,TYPE *f2,TYPE *f3, int N)
{
  int i;
  for(i=0;i<N;i++)
    {
      TYPE p = P(u1[i],u2[i],u3[i]);
      f1[i]=u1[i];
      f2[i]=p+pow(u2[i],2)/u1[i];
      f3[i]=(u2[i]/u1[i])*(u3[i]+p);
    }
}
