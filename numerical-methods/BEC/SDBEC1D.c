#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const size_t t=20000; //max iterations
const size_t n=400; //x points
const double x0=-16;
const double xn=16;
const double h=(xn-x0)/n;
#define TOL_RES 1e-12
//n=subintervalos, el nº de puntos es npts=n+1

double minimum(size_t n, double *arr)
{
    double min=arr[0];

    for (size_t i=0; i<n; i++)
    {
        if (arr[i]<min)
        {
            min=arr[i];
        }

    }

    return min;
}

void show_array(const double* array, size_t n)

{

    printf("[");

    for (int i = 0; i < n; i++) 
    {
        printf("%.7f", array[i]);
        
        if (i < n - 1) printf(", ");
    }

    printf("]\n");
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

double num_int_trap(size_t n, double a, double b, double *f) //metodo integral numérica por trapecios
{
    double sumf=0;

    double step=(b-a)/n;

    for (size_t i=1; i<n;i++)
    {
        sumf=sumf+f[i];
    }

    return step/2*(f[0]+2*sumf+f[n]);
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

double Energy(size_t n, double *xarr, double *V, double *phiold, double alpha, double *z) //Calculo de la energía para el estado fundamental dada la GPE a partir de una serie de valores para phi
{
    double GradPhi[n+1], fint[n+1], phi[n+1];

    for (size_t i=0; i<n+1; i++) phi[i]=phiold[i]-alpha*z[i];

    double nrm=num_norm_trap(n, x0, xn, phi);

    phi[0]=0; phi[n]=0;

    for (size_t i=0; i<n+1; i++) phi[i]/=nrm;

    num_grad(n,phi, GradPhi);

    for (size_t i=0; i<n+1; i++) fint[i]=0.5*GradPhi[i]*GradPhi[i]+V[i]*phi[i]*phi[i]+0.5*beta()*phi[i]*phi[i]*phi[i]*phi[i];

    double E=num_int_trap(n, x0, xn, fint);

    return E;
}

int main()
{
    double xarr[n+1], phi[n+1], phiOld[n+1], V[n+1], GradPhi[n+1], fint[n+1], LapPhi[n+1], GradPhi3[n+1], GradPhi2[n+1], GradPhi0[n+1], gradE[n+1], g3[n+1], g2[n+1], g0[n+1];

    double *E=malloc((t+1)*sizeof(double)), *E3=malloc((t+1)*sizeof(double)), *E2=malloc((t+1)*sizeof(double)), *E0=malloc((t+1)*sizeof(double));

    linspace(x0,xn,n,xarr);

    V_I(n, xarr,V);

    phi0(n,xarr, phi);

    size_t i;

    for (i=0; i<t+1; i++)
    {
        double alpha1=0;

        num_lap(n, phi, LapPhi);

        for(size_t i=0; i<n+1; i++)
        {
            gradE[i]=(-0.5*LapPhi[i]+V[i]*phi[i]+beta()*phi[i]*phi[i]*phi[i]);
        }

        double NormGrad=num_norm_trap(n, x0, xn, gradE);
        
        if (NormGrad==0)
        {
            printf("Zero gradient");

            break;
        }

        for (size_t i=0; i<n+1; i++) gradE[i]/=NormGrad; //z=z/z0

        E[i]=Energy(n, xarr, V, phi, alpha1, gradE);

        double alpha3=1;

        E3[i]=Energy(n, xarr, V, phi, alpha3, gradE);

        while (E3[i] >= E[i]) //Bucle para reducir E3 hasta que sea menor que E, o si no es posible, hasta que alpha 3 sea mas pequeño que la tolerancia
        {
            alpha3 *= 0.5;

            E3[i]=Energy(n, xarr, V, phi, alpha3, gradE);

            if (alpha3 < TOL_RES*0.5) 
            {
                printf("Sin mejora posible en it=%zu (res=%.3e) E %.15f\n", i, NormGrad, E3[i]);
            
                return 0;
            }

        }

        double alpha2=alpha3*0.5;

        E2[i]=Energy(n, xarr, V, phi, alpha2, gradE);

        //Interpolacion del polinomio para encontrar otra posible energía mediante alpha0

        double h1=(E2[i]-E[i])/alpha2;

        double h2=(E3[i]-E2[i])/(alpha3-alpha2);

        double h3=(h2-h1)/alpha3;

        double alpha0=0.5*(alpha2-h1/h3);

        E0[i]=Energy(n, xarr, V, phi, alpha0, gradE);

        //Buscamos el mínimo entre las posibles energías calculadas

        double alpha[4]={alpha0, alpha1, alpha2, alpha3};

        double tryE[4]={E0[i],E[i],E2[i],E3[i]};

        double min=minimum(4, tryE);

        if (fabs(E[i] - min) < TOL_RES) break;    // criterio para parar: |E - E_min| < Tol

        double alphaN;

        for (size_t i=0; i<4; i++)
        {
            alphaN=alpha[i];

            if (min==tryE[i])
            {
                break;
            }
        }

        for (size_t i=0; i<n+1; i++)
        {
            phiOld[i]=phi[i];
        }

        for (size_t i=0; i<n+1; i++)
        {
            phi[i]=phiOld[i]-alphaN*gradE[i];
        }

        double nrm=num_norm_trap(n, x0, xn, phi);

        phi[0]=0; phi[n]=0;

        for (size_t i=0; i<n+1; i++) phi[i]/=nrm;
    }
 
    printf("iteraciones=%zu  E_final=%.10f\n", i, E[i>t?t:i]);

    return 0;
}