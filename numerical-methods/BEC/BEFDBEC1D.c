#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const size_t t=20000; //max iterations
const size_t n=400; //x points
const double x0=-16;
const double xn=16;
const double h=(xn-x0)/n;
const double tao=0.1;
#define TOL_RES 1e-12

double maximum(size_t n, double *arr)
{
    double max=arr[0];

    for (size_t i=0; i<n; i++)
    {
        if (arr[i]>max)
        {
            max=arr[i];
        }

    }

    return max;
}

void linspace(double x0, double xn, size_t n, double *arr) //genera un array de n puntos equiespaciados en el intervalo x0 a xn
{
    double step = (xn-x0)/(n);
    
    for (int i = 0; i <= n; i++) 
    {
        arr[i] = x0 + i*step;
    }
}

void V_I(size_t n, double *xarr, double *V) //potencial usado para la simulación
{
    double gammax=1;

    for (size_t i=0; i<n+1; i++)
    {
        V[i]=0.5*gammax*gammax*xarr[i]*xarr[i];
    }
}

double beta() //cte beta
{
    double m=1.44e-25;

    double wx=20*M_PI,wy=20*M_PI,wz=20*M_PI;

    double hplanck=1.05e-34;

    double gammay=wy/wx;

    double gammaz=wz/wx;

    double xs=sqrt(hplanck/(m*wx));

    double N=133592.87;

    double as=5.1e-9;

    double kappa=(4*M_PI*as*N)/xs;

    return kappa*sqrt(gammay*gammaz)/(2*M_PI);
}

void phi0(size_t n, double *x, double *phi0) //aproximación inicial utilizada
{
    for (size_t i=0; i<n+1; i++)
    {
        phi0[i]=pow(M_PI,-0.25)*exp(-(x[i]*x[i])/2);
    }
}

void num_grad(size_t n, double *f, double *grad) //gradiente (derivada) numérica
{  
    grad[0]=0;

    grad[n]=0;

    for (size_t i=1; i<n; i++)
    {
        grad[i]=(f[i+1]-f[i])/(h);
    }
}

void num_lap(size_t n, double *f, double *lap) //laplaciano (derivada segunda) numerico
{
    lap[0]=0;

    lap[n]=0;

    for (size_t i=1; i<n; i++)
    {
        lap[i]=(f[i+1]-2*f[i]+f[i-1])/(h*h);
    }
}

double num_norm_trap(size_t n, double a, double b, double *f) //metodo para calcular la norma mediante la integral por trapecios
{
    double sumf=0;

    double step=(b-a)/n;

    for (size_t i=1; i<n;i++)
    {
        sumf=sumf+f[i]*f[i];
    }

    return sqrt(step/2*(f[0]*f[0]+2*sumf+f[n]*f[n]));
}

int  main()
{
    //x_(i+1)=(rho_i-x_i)/gamma_i: gamma_i=c_i/(b_i-a_i*gamma_(i-1)): rho_i=(r_i-a_i*rho_(i-1))/(b_i-a_i*gamma_(i-1))
    
    double gamma[n+1], rho[n+1], phi[n+1], xarr[n+1], b[n+1], V[n+1], phiold[n+1], deltaphi[n+1];

    double a=-tao/(2*h*h);

    double B=beta();

    linspace(x0,xn,n,xarr);

    V_I(n, xarr,V);

    phi0(n,xarr, phi);

    phi[0]=0; phi[n]=0;

    gamma[0]=a/(1+tao/(h*h)+tao*(V[0]+B*phi[0]*phi[0]));

    rho[0]=phi[0]/(1+tao/(h*h)+tao*(V[0]+B*phi[0]*phi[0]));

    size_t i;

    for (i=0; i<t+1; i++)
    {
        for (size_t j=0; j<n+1; j++) phiold[j]=phi[j];

        gamma[0]=0; rho[0]=0;

        for (size_t j=1; j<n; j++)
        {
            b[j]=1+tao/(h*h)+tao*(V[j]+B*phi[j]*phi[j]);

            gamma[j]=a/(b[j]-a*gamma[j-1]);

            rho[j]=(phi[j]-a*rho[j-1])/(b[j]-a*gamma[j-1]);
        }

        phi[n]=0;

        phi[n-1]=rho[n-1];

        for (size_t j=n-2; j>=1; j--) phi[j]=rho[j]-gamma[j]*phi[j+1];
 
        phi[0]=0;
 
        double nrm=num_norm_trap(n, x0, xn, phi);
 
        for (size_t k=0; k<n+1; k++) phi[k]/=nrm;

        for (size_t k=0; k<n+1; k++) deltaphi[k]=fabs(phi[k]-phiold[k]);
        
        double max=maximum(n, deltaphi);

            if (max < TOL_RES)
            {
                double E, gradPhi[n+1];

                num_grad(n, phi, gradPhi);

                for (size_t k=0; k<n+1; k++) E+=0.5*gradPhi[k]*gradPhi[k]+V[k]*phi[k]*phi[k]+B/2*phi[k]*phi[k]*phi[k]*phi[k];

                printf("iteraciones=%zu  E_final=%.10f\n", i, h*E);

                return 0;
            }
    }

    double E, gradPhi[n+1];

    num_grad(n, phi, gradPhi);

    for (size_t j=0; j<n+1; j++) E+=0.5*gradPhi[j]*gradPhi[j]+V[j]*phi[j]*phi[j]+B/2*phi[j]*phi[j]*phi[j]*phi[j];

    printf("Maximum iterations exceded %zu, %.10f", t, E);

    return 0;
}