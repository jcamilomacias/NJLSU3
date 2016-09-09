#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multifit.h>
#include "NJLSU3.h"
#define PI 3.1416
#define m_ud 0.0055
#define m_s 0.1407
#define A 0.6023
#define Gs 1.835/(A*A)
#define K 12.36/pow(A,5)
#define Nc 3
#define Nf 3

int
   print_state (size_t iter, gsl_multiroot_fsolver * s)
   {
   printf ("iter = %3zu x = % .8f % .8f % .8f"
   "f(x) = % .3e % .3e % .3e\n",
   iter,
   gsl_vector_get (s->x, 0),
   gsl_vector_get (s->x, 1),
   gsl_vector_get (s->x, 2),
   gsl_vector_get (s->f, 0),
	   gsl_vector_get (s->f, 1),
	   gsl_vector_get (s->f, 2));
   }
struct rparams
{
  double T;
  double mu_u;
  double mu_d;
  double mu_s;
        
};
   
int main(void)
{
  double x_init[2],phi_u,phi_d,phi_s,parametros_u[10],parametros_s[10],resultI1_u,resultI1_s,
  parametros_d[10],resultI2_u,resultI1_d,resultI2_d,resultI2_s,errorI2_d,errorI1_u,errorI1_d,errorI1_s,
  errorI2_u,errorI2_s,F0,F,phi_u0,phi_s0;
 int l_u,l_s,j,l_d; //POTENCIAL QUIMICO
  double P;
  
  FILE *Mvsmu,*Pgraf;
  Mvsmu=fopen("Mvsmu_NJLSU3.dat","w");
  Pgraf=fopen("Presure_NJLSU3.dat","w");
  
  gsl_integration_workspace * w
  =gsl_integration_workspace_alloc (100000);
  
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  
  int status;
  size_t iter=0;
  
  const size_t n=3;
  gsl_vector *x=gsl_vector_alloc (n);
  
  x_init[0]=0.370;
   x_init[1]=0.360;
  x_init[2]=0.550;
  
  gsl_vector_set (x,0,x_init[0]);
  gsl_vector_set (x,1,x_init[1]);
  gsl_vector_set (x,2,x_init[2]);
  struct rparams p={};
  gsl_multiroot_function f = {&gap0, n, &p};
 
  T=gsl_multiroot_fsolver_dnewton;
  s=gsl_multiroot_fsolver_alloc (T,3);
  gsl_multiroot_fsolver_set (s,&f,x);
  
  do
  {
    iter++;
    status=gsl_multiroot_fsolver_iterate (s);
    
    if (status)
      break;
    
    status=gsl_multiroot_test_residual (s->f, 1e-12);
  }
  while (status == GSL_CONTINUE && iter < 1000);
  printf("% lf % lf\n",gsl_vector_get(s->x,0),gsl_vector_get(s->x,1));
  iter=0;
//x_init[0]=0.400;
//x_init[1]=0.600;
  x_init[0]=gsl_vector_get (s->x,0);
  x_init[1]=gsl_vector_get (s->x,1);
  x_init[2]=gsl_vector_get (s->x,2);
  
  phi_u0=-2*Nc*x_init[1]*I20(x_init[1]);
  phi_s0=-2*Nc*x_init[2]*I20(x_init[2]);
  
  F0=-2*Nc*2*I10(x_init[0])-2*Nc*I10(x_init[2])+2*Gs*(2*pow(phi_u0,2)+pow(phi_s0,2))-
  4*K*pow(phi_u0,2)*phi_s0;
  //F0=-2*Nc*2*I10(x_init[0])+2*Gs*(2*pow(phi_u0,2));
  gsl_vector_set (x,0,x_init[0]);
  gsl_vector_set (x,1,x_init[1]);
  gsl_vector_set (x,2,x_init[2]);
  printf("%lf\n",F0);
  
  for(j=1;j<=400;j++)
  {
    for(l_u=0;l_u<=50;l_u++)
    {
      for(l_d=0;l_d<=50;l_d++)
      {
      for(l_s=0;l_s<=5;l_s++)
	{
	  double j_T,mu_u,mu_d,mu_s;
	  j_T=j/1000.0;
	  mu_u=l_u/1000.0;
	  mu_d=l_d/1000.0;
	  mu_s=l_s/1000.0;
	struct rparams p = {j_T,mu_u,mu_d,mu_s};
	gsl_multiroot_function f={&gap,n,&p};
    
	T=gsl_multiroot_fsolver_dnewton;
	s= gsl_multiroot_fsolver_alloc (T,3);
	gsl_multiroot_fsolver_set (s,&f,x);
    
	do
	{
	  print_state (iter,s);
	  iter++;
	  status = gsl_multiroot_fsolver_iterate (s);
      
	  if (status)
	  break;
      
	  status = gsl_multiroot_test_residual (s->f, 1e-12);
	}
    while (status == GSL_CONTINUE && iter < 1000);
    iter=0;
    
    
  
  fprintf(Mvsmu,"%lf %lf %lf %lf\n",j_T,gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),gsl_vector_get(s->x,2));
   x_init[0]=gsl_vector_get (s->x,0);
   x_init[1]=gsl_vector_get (s->x,1);
  x_init[2]=gsl_vector_get (s->x,2);
   
  parametros_u[0]=x_init[0];
  parametros_u[1]=j_T;
  parametros_u[2]=mu_u;

  parametros_d[0]=x_init[1];
  parametros_d[1]=j_T;
  parametros_d[2]=mu_d;
  
  parametros_s[0]=x_init[2];
  parametros_s[1]=j_T;
  parametros_s[2]=mu_s;
  
  gsl_function Q_u; 
      Q_u.function = &I1;
      Q_u.params = &parametros_u;
      gsl_integration_qagiu (&Q_u, 0, 0, 1e-12, 100000,
			 w, &resultI1_u, &errorI1_u); 
      gsl_function Q_d; 
      Q_d.function = &I1;
      Q_d.params = &parametros_d;
      gsl_integration_qagiu (&Q_d, 0, 0, 1e-12, 100000,
			 w, &resultI1_d, &errorI1_d); 
      
      gsl_function Q_s; 
      Q_s.function = &I1;
      Q_s.params = &parametros_s;
      gsl_integration_qagiu (&Q_s, 0, 0, 1e-12, 100000,
			 w, &resultI1_s, &errorI1_s);
     gsl_function R_u; 
      R_u.function = &I2;
      R_u.params = &parametros_u;
      gsl_integration_qagiu (&R_u, 0, 0, 1e-12, 100000,
			 w, &resultI2_u, &errorI2_u);

	gsl_function R_d; 
      R_d.function = &I2;
      R_d.params = &parametros_d;
      gsl_integration_qagiu (&R_d, 0, 0, 1e-12, 100000,
			 w, &resultI2_d, &errorI2_d);
      
      gsl_function R_s; 
      R_s.function = &I2;
      R_s.params = &parametros_s;
      gsl_integration_qagiu (&R_s, 0, 0, 1e-12, 100000,
			 w, &resultI2_s, &errorI2_s);
      
      phi_u=-2*Nc*x_init[0]*(I20(x_init[0])+resultI2_u);
      phi_d=-2*Nc*x_init[1]*(I20(x_init[1])+resultI2_d);
     // printf("%lf\n",phi_u);
      phi_s=-2*Nc*x_init[2]*(I20(x_init[2])+resultI2_s);
 F=-2*Nc*(I10(x_init[0])+resultI1_u)-2*Nc*(I10(x_init[1])+resultI1_d)-2*Nc*(I10(x_init[2])+resultI1_s)
 +2*Gs*(pow(phi_u,2)+pow(phi_d,2)+pow(phi_s,2))-4*K*phi_u*phi_d*phi_s;
      // F=-2*Nc*2*(I10(x_init[0])+resultI1_u)-2*Nc*(I10(x_init[2])+resultI1_s)
//+2*Gs*(2*pow(phi_u,2)+pow(phi_s,2))-4*K*pow(phi_u,2)*phi_s;
 //F=-2*Nc*2*(I10(x_init[0])+resultI1_u)
 //+2*Gs*(2*pow(phi_u,2));
// P[400*4001*6*(j-1)+401*6*l_u+6*l_d+l_s]=(float)(F0-F)/pow(j,4);
 P=(F0-F)/pow(j_T,4);
 fprintf(Pgraf,"%lf %lf %lf % lf % lf\n",j_T,mu_u,mu_d,mu_s,P);
 
  gsl_vector_set (x,0,x_init[0]);
  gsl_vector_set (x,1,x_init[1]);
  gsl_vector_set (x,2,x_init[2]);
	}
    }
  }
   printf("T=%lf\n",j);
  }
   gsl_integration_workspace_free (w);   
  gsl_multiroot_fsolver_free (s);
   gsl_vector_free (x);
   fclose(Mvsmu);
   fclose(Pgraf);
   return 0;
}