from sympy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from sympy.plotting import plot
from matplotlib.widgets import TextBox, Slider
mpl.use('TkAgg')
from typing import Tuple
from tabulate import tabulate
def ajuste_linear(x:NDArray[np.float64],y:NDArray[np.float64],errorx:NDArray[np.float64] | float=0,errory:NDArray[np.float64] | float=0,xlabel:str='',ylabel:str='',titulo:str='',color:str='orange',colorerr:str='black',colordot:str='red'):
    """
    Realiza un ajuste lineal de los datos x e y introducidos en la función, y devuelve la pendiente, la ordenada en origen y el coeficiente de correlación al cuadrado.
    """
    n=len(x)
    xm=np.mean(x)
    ym=np.mean(y)
    sx=np.sum(x)
    sy=np.sum(y)
    sxy=np.sum(x*y)
    sx2=np.sum(x**2)
    sy2=np.sum(y**2)
    m=(n*sxy-sx*sy)/(n*sx2-sx**2)
    b=(sx2*sy-sx*sxy)/(n*sx2-sx**2)
    xi=symbols('xi')
    f=b+m*xi
    fx=lambdify(xi,f)
    fi=fx(x)
    numerador=n*sxy-sx*sy
    raiz1=np.sqrt(n*sx2-sx**2)
    raiz2=np.sqrt(n*sy2-sy**2)
    r=numerador/(raiz1*raiz2)
    r2=r**2
    r2_porcentaje=np.around(r2*100,2)
    sigma=(sum((y-m*x-b)**2)/(len(x)-2))**(1/2)
    errorm=(((len(x)*sigma**2)/(len(x)*sum(x**2)-(sum(x))**2)))**(1/2)
    errorn=((sigma**2*sum(x**2))/(len(x)*sum(x**2)-(sum(x))**2))**(1/2)
    m1='m='+str(m)+'+/-'+str(errorm)
    b1='n='+str(b)+'+/-'+str(errorn)
    plt.plot(x,y,'o',label='Datos')
    plt.plot(x,fi,color=color,label='Ajuste')
    plt.errorbar(x,y,color=colorerr,yerr=errory,fmt='.')
    plt.errorbar(x,y,color=colorerr,xerr=errorx,fmt='.')
    for i in range(0,n,1):
        y0=np.min([y[i],fi[i]])
        y1=np.max([y[i],fi[i]])
        plt.vlines(x[i],y0,y1, color=colordot,
            linestyle ='dotted')
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(titulo)
    plt.show
    f1='f='+str(f)
    r21='r^2='+str(r2)
    m1val=m
    return f1,m1,b1,r21,m1val

def properrores(variables:list[Symbol],f:Symbol,errores:list[float],valores1:NDArray[np.float64]):
    """
    Es necesario definir de manera simbólica antes de realizar la propagación de errores.
    En los valores se ha de utilizar un array bidimensional (una matriz).
    Realiza la propagación de errores de una función dadas unas variables y unos valores con sus respectivos errores, para devolver la funcion de error y los valores de estos.
    """
    ecuacion=0
    errorindirecto1=[]
    valores=valores1.T
    for i in range(len(variables)):
        ecuacion+=abs(diff(f,variables[i]))*errores[i]
    print(ecuacion)
    eqlatex=latex(ecuacion)
    print(eqlatex)
    for j in range(len(valores)):
        errorindirecto2=ecuacion.subs(dict(zip(variables,valores[j])))
        errorindirecto1.append(errorindirecto2)
    return np.array(errorindirecto1)
                #mA    #mV
Tabla1=np.array([[199.2,36.6],
                [194.5,35.6],
                [189.8,34.7],
                [184.7,33.8],
                [180.1,32.9],
                [174.5,31.9],
                [169.9,31.0],
                [165.1,30.2],
                [159.9,29.2],
                [154.5,28.3],
                [149.9,27.4],
                [144.6,26.4],
                [139.8,25.6],
                [135.4,24.7],
                [129.0,23.7],
                [124.5,22.7],
                [119.8,21.9],
                [115.0,21.0],
                [110.1,20.1],
                [105.6,19.3],
                [100.1,18.3]]
                )
print(tabulate(Tabla1,['Intensidad','Voltaje'],tablefmt="latex"))
a1=ajuste_linear(Tabla1.T[0],Tabla1.T[1],0.1,0.1,xlabel='I[mA]',ylabel='V[mV]',titulo='I frente a V')
m1=a1[4]
plt.show()
print(a1)
print(m1)
Tabla3=np.array([[0.28,0.05],
                 [0.38,0.07],
                 [0.47,0.09],
                 [0.57,0.11],
                 [0.67,0.13],
                 [0.77,0.15],
                 [0.82,0.16],
                 [0.88,0.18],
                 [0.99,0.21],
                 [1.08,0.23],
                 [1.20,0.27],
                 [1.31,0.31],
                 [1.40,0.35],
                 [1.47,0.39],
                 [1.58,0.47]])
a3=ajuste_linear(Tabla3.T[0],Tabla3.T[1],0.01,0.01)
m4=a3[4]
print(a3)
print(m4)
plt.show()
                #U[V]#I[A]#U_therm[mV] #Subiendo en pasos de 0.5V empezando en 1V
Tabla2=np.array([[1,1.99,0.23],
                 [1.43,2.32,0.34],
                 [1.89,2.63,0.58],
                 [2.43,2.96,0.93],
                 [2.93,3.25,1.33],
                 [3.4,3.51,1.72],
                 [3.86,3.74,2.15],
                 [4.48,4.00,2.70],
                 [4.81,4.20,3.18],
                 [5.32,4.43,3.84],
                 [5.87,4.67,4.50]])
print(tabulate(Tabla2,['Intensidad','Voltaje'],tablefmt="latex"))
a2=ajuste_linear(Tabla2.T[0],Tabla2.T[1],0.01,0.01)
m2=a2[4]
plt.show()
print(a2)
print(m2)
tR=296
alpha=4.82*10**-3
beta=6.76*10**-7
RtR=m4
Rt=Tabla2.T[0]/Tabla2.T[1]
R0=RtR/(1+alpha*tR+beta*tR**2)
T=273+(1/(2*beta))*(np.sqrt(alpha**2+4*beta*(Rt/R0-1))-alpha)
print(T)
a3=ajuste_linear(np.log(T),np.log(Tabla2.T[2]))
m3=a3[4]
print(a3)
print(m3)
plt.show()
#Calcular dato con Tabla 1 y Tabla 3
#V1=5V
               #D1[cm] U_therm1[mV] D2[cm] U_therm2[mV]
Tabla4=np.array([[20,5.84,20,7.14],
                 [25,3.30,25,4.37],
                 [30,1.87,30,2.53],
                 [35,1.16,35,1.64],
                 [40,0.78,40,0.98],
                 ])
