#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TYPE double

void inicializar(TYPE *p,TYPE *rho, TYPE *v,TYPE *x,int N);

int main()
{
  TYPE delta_x = 0.01;
  TYPE delta_t=0.1;
  TYPE l=1.0;
  
  int N=l/delta_x;

  TYPE* P= malloc(N*sizeof(TYPE));
  TYPE* rho= malloc(N*sizeof(TYPE));
  TYPE* v= malloc(N*sizeof(TYPE));
  TYPE* x= malloc(N*sizeof(TYPE));

  inicializar(p,rho,v,x,N,delta_x);
  
}


void inicializar(TYPE *p,TYPE *rho, TYPE *v,TYPE *x,int N,TYPE dx )
{
  int Naux=0.5/dx 
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
