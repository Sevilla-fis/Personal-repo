from typing import Callable
from sympy import *
import matplotlib as mpl
from matplotlib import animation
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from statistics import mean, pstdev
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from sympy.plotting import plot
from matplotlib.widgets import TextBox, Slider
mpl.use('TkAgg')
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
#DIA 1
#PARTE 1. MASAS IGUALES. SIN PESAS
l=0.1
t_1i=np.array([0.159,0.147,0.119,0.084,0.172,0.112])
t_1f=np.array([1.432,1.364,0.851,0.612,0.843,0.562])
t_2=np.array([0.237,0.207,0.176,0.162,0.261,0.192])
vi_1=l/t_1i
vf_2=l/t_2
rep=ajuste_linear(vi_1,vf_2)
plt.show()
print(rep)
#Masa del carrito sin pesas
m=0.19051#Kg
#PARTE 2. LE AÑADIMOS IMANES
m=0.19718
t_1i_1=np.array([0.199,0.233,0.305,0.219,0.200,0.296,0.224,0.254,0.243,0.244])
t_2f_1=np.array([0.311,0.336,0.383,0.316,0.302,0.365,0.322,0.338,0.326,0.335])
vi_1_1=l/t_1i_1
vf_2_1=l/t_2f_1
rep_1=ajuste_linear(vi_1_1,vf_2_1)
print(rep_1)
print(vi_1_1,vf_2_1)
plt.show()
#PARTE 3. LE AÑADIMOS MASA(100g) AL CARRO 2 ADEMAS
t_1i_2=np.array([0.315,0.336,0.292,0.233,0.244,0.223,0.257,0.190,0.258,0.218])
t_2f_2=np.array([0.509,0.551,0.516,0.488,0.507,0.492,0.482,0.385,0.488,0.470])
vi_1_2=l/t_1i_2
vf_2_2=(l/t_2f_2)
rep_2=ajuste_linear(vi_1_2,vf_2_2)
plt.show()
print(rep_2)
#DIA 2
l=0.1
t_1i_1=np.array([0.214,0.233,0.209,0.164,0.195,0.186,0.167,0.142,0.154,0.144,0.179])
t_2f_1=np.array([0.221,0.243,0.218,0.171,0.203,0.195,0.174,0.150,0.161,0.150,0.187])
v1i_1=l/t_1i_1
v2f_1=l/t_2f_1
print(ajuste_linear(v1i_1,v2f_1))
plt.show()
t_1i_2=np.array([0.166,0.159,0.180,0.152,0.179,0.180,0.162,0.187,0.169,0.155,0.169])
t_2f_2=np.array([0.431,0.392,0.434,0.378,0.429,0.446,0.410,0.464,0.428,0.393,0.422])
v1i_2=l/t_1i_2
v2f_2=l/t_2f_2
print(ajuste_linear(v2f_2,v1i_2))
plt.show()