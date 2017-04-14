#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TYPE float
#define dx 0.001


void inicializar(TYPE p[], TYPE rho[], TYPE v[], TYPE x[], int N);
void inicialUF(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], TYPE rho[], TYPE v[], TYPE p[], int N);
TYPE e(TYPE rho,TYPE p);
TYPE E(TYPE v,TYPE rho,TYPE e);
void LaxWendroff(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], int N);
TYPE P(TYPE u1,TYPE u2,TYPE u3);
void actualizaF(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], int N);
TYPE Umax(TYPE u1[], TYPE u2[], TYPE u3[], int N);

int main()
{
  int l = 1;
  int N = l/dx;
  
  TYPE p[N];
  TYPE rho[N];
  TYPE v[N];
  TYPE x[N];
  
  inicializar(p,rho,v,x,N);
  
  TYPE u1[N];
  TYPE u2[N];
  TYPE u3[N];
  TYPE f1[N];
  TYPE f2[N];
  TYPE f3[N];

  TYPE max;

  inicialUF(u1,u2,u3,f1,f2,f3,rho,v,p,N);
  int t=0;
  while(abs((u2[(int)(0.9/dx)]/u1[(int)(0.9/dx)])-(u2[(int)(0.9/dx+1)]/u1[(int)(0.9/dx+1)])) < 1E-5)   
    {
      int i;
      //printf("U y F antes de Lax");
      //for(i=0;i<N;i++)
      //	 {
      //	   printf("%f %f %f %f %f %f\n",u1[i],u2[i],u3[i],f1[i],f2[i],f3[i]);
      //	 }
      LaxWendroff(u1,u2,u3,f1,f2,f3,N);
      //printf("U y F despues de Lax");
      //for(i=0;i<N;i++)
      //	 {
      //	   printf("%f %f %f %f %f %f\n",u1[i],u2[i],u3[i],f1[i],f2[i],f3[i]);
      //	 }
      actualizaF(u1,u2,u3,f1,f2,f3,N);
      //printf("U y F despues de actualizarF");
      //for(i=0;i<N;i++)
      //	 {
      //	   printf("%f %f %f %f %f %f\n",u1[i],u2[i],u3[i],f1[i],f2[i],f3[i]);
      //	 }
      
      //for(i=0;i<N;i++)
      //	 {
      //	   printf("%f %f %f %f\n",x[i],u1[i],u2[i],u3[i]);
      //	 }
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
void inicializar(TYPE p[],TYPE rho[], TYPE v[],TYPE x[], int N)
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
	  rho[i]=0.125;
	}
      x[i]=i*dx;
      v[i]=0;
    }
}

/* Inicializa U y F */
void inicialUF(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], TYPE rho[], TYPE v[], TYPE p[], int N)
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

/* realiza el metodo lax-wendroff */
void LaxWendroff(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], int N)
{
  TYPE *tempu1= malloc(N*sizeof(TYPE));
  TYPE *tempu2= malloc(N*sizeof(TYPE));
  TYPE *tempu3= malloc(N*sizeof(TYPE));

  TYPE U1H1;
  TYPE U2H1;
  TYPE U3H1;
  TYPE U1H2;
  TYPE U2H2;
  TYPE U3H2;

  TYPE F1H1;
  TYPE F2H1;
  TYPE F3H1;
  TYPE F1H2;
  TYPE F2H2;
  TYPE F3H2;
  
  tempu1[0]=u1[0];
  tempu2[0]=u2[0];
  tempu3[0]=u3[0];

  tempu1[(int)(N-1)]=u1[(int)(N-1)];
  tempu2[(int)(N-1)]=u2[(int)(N-1)];
  tempu3[(int)(N-1)]=u3[(int)(N-1)];

  /* Inicial step */
  
  TYPE dt=0.01;
  int i;
  for (i=1;i<N-1;i++)
    {  
      U1H1=0.5*(u1[i+1]+u1[i])-0.5*(dt/dx)*(f1[i+1]-f1[i]);
      U2H1=0.5*(u2[i+1]+u2[i])-0.5*(dt/dx)*(f2[i+1]-f2[i]);
      U3H1=0.5*(u3[i+1]+u3[i])-0.5*(dt/dx)*(f3[i+1]-f3[i]);
      U1H2=0.5*(u1[i]+u1[i-1])-0.5*(dt/dx)*(f1[i]-f1[i-1]);
      U2H2=0.5*(u2[i]+u2[i-1])-0.5*(dt/dx)*(f2[i]-f2[i-1]);
      U3H2=0.5*(u3[i]+u3[i-1])-0.5*(dt/dx)*(f3[i]-f3[i-1]);

      F1H1=U2H1;
      F2H1=P(U1H1,U2H1,U3H1)+U2H1*U2H1/U1H1;
      F3H1=(U2H1/U1H1)*(U3H1+P(U1H1,U2H1,U3H1));
      F1H2=U2H2;
      F2H2=P(U1H2,U2H2,U3H2)+U2H2*U2H2/U1H2;
      F3H2=(U2H2/U1H2)*(U3H2+P(U1H2,U2H2,U3H2));

      tempu1[i]=u1[i]-(dt/dx)*(F1H1-F1H2);
      tempu2[i]=u2[i]-(dt/dx)*(F2H1-F2H2);;
      tempu3[i]=u3[i]-(dt/dx)*(F3H1-F3H2);;
    }
  
  int r;
  for(r=0;r<N;r++)
    {
       u1[r]=tempu1[r];
       u2[r]=tempu2[r];
       u3[r]=tempu3[r];
    }
  /* Falta mirar la condicion dt/dx umax<1 */
}


/* Halla la presion */

TYPE P(TYPE u1,TYPE u2,TYPE u3)
{
  return 0.4*(u3-0.5*u2*u2/(u1));
}

/* Actualiza el vector F */

void actualizaF(TYPE u1[],TYPE u2[],TYPE u3[],TYPE f1[],TYPE f2[],TYPE f3[], int N)
{
  int i;
  for(i=0;i<N;i++)
    {
      TYPE p = P(u1[i],u2[i],u3[i]);
      f1[i]=u2[i];
      f2[i]=p+u2[i]*u2[i]/u1[i];
      f3[i]=(u2[i]/u1[i])*(u3[i]+p);
    }
}


