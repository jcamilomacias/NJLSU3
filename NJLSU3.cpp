#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multifit.h>
#define PI 3.1416
#define m_ud 0.0055
#define m_s 0.1407
#define A 0.6023
#define Gs 1.835/(A*A)
#define K 12.36/pow(A,5)
#define Nc 3
#define Nf 3



double fplus(double x, void * parametros){
  double *p,E,mu,M,T;
  
  p=(double *) parametros;
  
  M=p[0];
  T=p[1];
  mu=p[2];
  
  E=sqrt(x*x+M*M);
  
  double fplus= exp(-(E-mu)/T)/(1+exp(-(E-mu)/T));
  return fplus;
}

double fminus(double x, void * parametros){
  double *p,E,mu,M,T;
  
  p=(double *) parametros;
  
  M=p[0];
  T=p[1];
  mu=p[2];
  
  E=sqrt(x*x+M*M);
  
  double fminus=exp(-(E+mu)/T)/(1+exp(-(E+mu)/T));
  return fminus;
}
double I10 (double x){
  double raiz;
  raiz=sqrt(A*A+pow(x,2));
  double I10=(-1/(32*PI*PI))*(pow(x,4)*log(pow(A+raiz,2)/pow(x,2))-2*raiz*(2*
    pow(A,3)+A*pow(x,2)));
  return I10;
}
double I1(double x, void * parametros ){
  double *p,M,T,E,arg1,mu,arg2;
   p = (double *) parametros;
   M=p[0];
  T=p[1];
  mu=p[2];
   E=sqrt(x*x+pow(M,2));
    arg1=1+exp(-(E+mu)/T);
   arg2=1+exp(-(E-mu)/T);
  
   double I1= 1/(2*PI*PI)*T*x*x*(log(arg1)+log(arg2));
   return I1;   
}

double I20(double M){
  
  double I20=(1/(4*PI*PI))*(A*sqrt(A*A+M*M)
  -0.5*M*M*log(pow(A+sqrt(A*A+M*M),2)/(M*M)));
  return I20;
}

double I2(double x, void * parametros){
  double *p,E,M;
  
  p=(double *) parametros;
  
  M=p[0];
  
  E=sqrt(x*x+M*M);
  
  double I2=-1/(2*PI*PI*E)*x*x*(fplus(x,p)+fminus(x,p));
  return I2;
}
struct rparams
{
  double T;
  double mu_u;
  double mu_d;
  double mu_s;
        
};

int gap0(const gsl_vector *x, void *params, gsl_vector * f)
{
  const double x0=gsl_vector_get (x,0);
  const double x1=gsl_vector_get (x,1);
  const double x2=gsl_vector_get (x,2);

  const double y0=x0-m_ud-4*Gs*2*Nc*x0*I20(x0)
  -2*K*4*Nc*Nc*x1*x2*I20(x1)*I20(x2);
  const double y1=x1-m_ud-4*Gs*2*Nc*x1*I20(x1)
  -2*K*4*Nc*Nc*x0*x2*I20(x0)*I20(x2);
  const double y2=x2-m_s-4*Gs*2*Nc*x2*I20(x2)-2*K*4*Nc*Nc*x0*x1*I20(x0)*I20(x1);

  
  gsl_vector_set (f,0,y0);
  gsl_vector_set (f,1,y1);
  gsl_vector_set (f,2,y2);
  
  return GSL_SUCCESS;
}

int gap (const gsl_vector * x, void *params, gsl_vector * f)
{
  double integralI2_u, error_u,integralI2_s,error_s,integralI2_d,error_d,
  param[10],param2[10],param3[10];
  gsl_integration_workspace * w
  =gsl_integration_workspace_alloc (100000);
  
  double T= ((struct rparams *) params)->T;
  double mu_u= ((struct rparams *) params)->mu_u;
  double mu_d= ((struct rparams *) params)->mu_d;
  double mu_s= ((struct rparams *) params)->mu_s;
  
  const double x0=gsl_vector_get (x,0);
  const double x1=gsl_vector_get (x,1);
  const double x2=gsl_vector_get (x,2);
  
  param[0]=x0;
  param[1]=T;
  param[2]=mu_u;
   
 param2[0]=x1;
  param2[1]=T;
  param2[2]=mu_d;
  
  param3[0]=x2;
  param3[1]=T;
  param3[2]=mu_s;
  
  gsl_set_error_handler_off();
  gsl_function B;
  B.function=&I2;
  B.params=&param;
  gsl_integration_qagiu(&B,0,0,1e-12, 10000,w,&integralI2_u,&error_u);
  
  gsl_function C;
  C.function=&I2;
  C.params=&param2;
  gsl_integration_qagiu (&C,0,0,1e-12, 10000,w,&integralI2_d,&error_d);
  
  gsl_function D;
  D.function=&I2;
  D.params=&param3;
  gsl_integration_qagiu (&D,0,0,1e-12, 10000,w,&integralI2_s,&error_s);
  
  const double y0=x0-m_ud-4*Gs*2*Nc*x0*(I20(x0)+integralI2_u)
  -2*K*4*Nc*Nc*x1*x2*(I20(x1)+integralI2_d)*(I20(x2)+integralI2_s);
 const double y1=x1-m_ud-4*Gs*2*Nc*x1*(I20(x1)+integralI2_d)
  -2*K*4*Nc*Nc*x0*x2*(I20(x0)+integralI2_u)*(I20(x2)+integralI2_s);
  const double y2=x2-m_s-4*Gs*2*Nc*x2*(I20(x2)+integralI2_s)
  -2*K*4*Nc*Nc*x0*x1*(I20(x0)+integralI2_u)*(I20(x1)+integralI2_d);

gsl_vector_set (f,0,y0);
gsl_vector_set (f,1,y1);
gsl_vector_set (f,2,y2);

gsl_integration_workspace_free (w);
return GSL_SUCCESS;
}
