from typing import Callable
from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from statistics import mean, pstdev
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from sympy.plotting import plot
import matplotlib as mpl
mpl.use('TkAgg')
def csvtopython(csvfile,dtype=None,delimiter=';'):

    data=np.genfromtxt(csvfile,delimiter=delimiter,unpack=True,dtype=dtype) #We obtain data with np.genfromtxt
    #It's important that data is separated by . in excel (0.23), so type is float
    #Delimiter is set to ; (Change if necessary)

    datos=[]

    for i in range(len(data)): #Remove first element of array (name of variable) so we obtain only data (type:float)

        dat=data[i][1:]

        datos.append(dat)

    return datos

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
    return f1,m1,b1,r21

def ajuste_magico(x:NDArray[np.float64],y:NDArray[np.float64],funcion:Callable[...,float] ,parametros:list[str],p0:NDArray[np.float64],xerror:NDArray[np.float64]|float=0,yerror:NDArray[np.float64]|float=0,xlabel:str='',ylabel:str='',titulo:str='',color:str='r-',colorerr:str='black'):
    """
    Toma una función como argumento, por lo tanto es necesario definirla antes de realizar el ajuste.
    Realiza un ajuste de cualquier función y devuelve los parámetros de la misma con su error correspondiente.
    """
    popt,pcov=curve_fit(funcion,x,y,p0)
    min_x=x[0]
    for i in x[1:]:
        min_x=min(min_x,i)
    max_x=x[-1]
    for i in x[1:]:
        max_x=max(max_x,i)
    x_ajuste=np.linspace(min_x,max_x)
    plt.figure()
    plt.plot(x,y,'o',label='Datos')
    plt.plot(x_ajuste,funcion(x_ajuste,*popt),color,label='Ajuste')
    plt.errorbar(x,y,color=colorerr,yerr=yerror,fmt='.')
    plt.errorbar(x,y,color=colorerr,xerr=xerror,fmt='.')
    plt.legend(loc='best')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(titulo)
    plt.tight_layout()
    pstd=np.sqrt(np.diag(pcov))
    nombres_de_param=parametros
    print('Parámetros:')
    for i, param in enumerate(popt):
        print('{:s}={:f}+-{:f}'.format(nombres_de_param[i],param,pstd[i]/2))
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
B1=198*10**-3
d=10**-3
I1=10**-3*np.array([-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])
x1=I1*B1/d
Uh1=10**-3*np.array([-25.5,-20.17,-16.36,-11.53,-6.66,-2.02,2.19,6.07,11.32,15.89,20.14,24.37,29.35])
print(ajuste_linear(x1,Uh1))
plt.show()
I2=-31*10**-3
B2=10**-3*np.array([9,20,30,42,51,63,73,80,91,102,110,119,131,141,150,162,171,182,192,201])
Uh2=10**-3*np.array([-0.24,-1.84,-3.16,-4.87,-6.24,-7.89,-9.27,-10.29,-11.87,-13.36,-14.46,-15.73,-17.33,-18.74,-19.90,-21.63,-22.78,-24.20,-25.48,-26.65])
x2=I2*B2/d
print(ajuste_linear(x2,Uh2))
plt.show()
Uh3=10**-3*np.array([5.1,5.18,5.15,5.05,5.07,5.26,5.42,5.54,5.64,5.74,5.65,5.32,4.87,4.07,3,1.87,0.47,-0.62,-1.44,-1.88,-2.07,-2.09])
T=273+np.array([28,33,38,43,48,53,58,63,68,73,78,83,88,93,98,103,108,113,118,123,128,133])
plt.plot(T,Uh3,'o')
plt.show()
Usample=np.array([0.9590,0.9689,0.9895,1.0190,1.0546,1.0853,1.1126,1.1493,1.1760,1.2095,1.2216,1.2312,1.2261,1.1986,1.15654,1.1388,1.0065,0.9198,0.8150,0.7235,0.6272,0.5410])
plt.plot(1/T,-np.log(Usample),'o')
plt.show()
