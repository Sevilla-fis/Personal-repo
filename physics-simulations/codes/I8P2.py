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
#Angulo Brewster
plt.figure()
r_p=np.array([0.350,0.106,0.049,0.020,0.010,0.003,0.001,0.000,0.001,0.002,0.003,0.004,0.004,0.005,0.005,0.005])
r_s=np.array([0.350,0.095,0.064,0.030,0.026,0.025,0.020,0.015,0.019,0.015,0.013,0.011,0.008,0.006,0.003,0.004])
Angulo=np.pi/180*np.array([90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15])
plt.figure()
plt.plot(Angulo,r_p, 'o', label='Datos')
plt.errorbar(Angulo, r_p, yerr=0.001, fmt='.')
plt.errorbar(Angulo, r_p, xerr=np.pi/180, fmt='.')
plt.legend(loc='best')
plt.xlabel('Angulo [rad]')
plt.ylabel('r_p [W]')
plt.title('Intensidad paralela reflejada frente al ángulo')
plt.tight_layout()
plt.show()
I=np.array([0.014,0,0.009,0.040,0.082,0.147,0.203,0.264,0.310,0.337,0.358,0.347,0.312,0.263,0.210,0.148,0.090,0.046,0.010])
theta=-10+np.array([-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90])
theta1=symbols("theta")
f=tan(theta1)
err=properrores([theta1],f,[1*np.pi/180],np.array([[55*np.pi/180]]))
print('Error n_2=', err)
I0=np.array([0.35])
def coseno2(x,I0):
    return I0*np.cos(x)**2
ajuste_magico(theta*np.pi/180,I,coseno2,["I0"],I0,xerror=1*np.pi/180,yerror=0.001,xlabel="Ángulo(rad)",ylabel="Intensidad [W]",titulo="Ley de Malus")
plt.show()
ILinealVertical=np.array([0.54,0.35,0,0.177,0.177,0.296,0.288])
ICircularDextrogira=np.array([0.440,0.241,0.021,0.074,0.201,0.440,0.140])
"np.array([Itotal,IV,IH,I45,I-45,ICD,ICI"
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1,1,1,projection='3d')
ax.set_xlim3d(0,100)
ax.set_ylim3d(-2,2)
ax.set_zlim3d(-2,2)
fig.subplots_adjust(bottom=0.15)
axslider1 = fig.add_axes([0.09,0.01,0.2,0.1])
axslider2 = fig.add_axes([0.4,0.01,0.2,0.1])
axslider3 = fig.add_axes([0.7,0.01,0.2,0.1])
sliderE0x = Slider(axslider1, "E0x",0,1)
sliderE0y = Slider(axslider2, "E0y",0,1)
sliderPhi = Slider(axslider3, "Phi",0,2*np.pi)
w=1
t = np.linspace(0,100,1000)
Ex=sliderE0x.val * np.cos(w * t)
Ey=sliderE0y.val * np.cos(w * t+sliderPhi.val)
line = ax.plot(t[0],Ex[0],lw=2)[0]
def animate(i):
    Ex = sliderE0x.val * np.cos(w * t)
    Ey = sliderE0y.val * np.cos(w * t + sliderPhi.val)
    x = t[:i]
    y = Ex[:i]
    z = Ey[:i]
    s0=sliderE0x.val ** 2 + sliderE0y.val ** 2
    s1=sliderE0x.val ** 2 - sliderE0y.val ** 2
    s2=2 * sliderE0x.val * sliderE0y.val * np.cos(sliderPhi.val)
    s3=2 * sliderE0x.val * sliderE0y.val * np.sin(sliderPhi.val)
    norm=np.sqrt(s0**2+s1**2+s2**2+s3**2)
    plt.legend(["S0={}".format(s0/norm),
                "S1={}".format(s1/norm),
                "S2={}".format(s2/norm),
                "S3={}".format(s3/norm )],bbox_to_anchor=(1.4, 9.5), ncol=1)
    line.set_data(x,y)
    line.set_3d_properties(z)
    return line
anim = animation.FuncAnimation(fig=fig,func=animate, frames=10000, interval=10)
plt.show()