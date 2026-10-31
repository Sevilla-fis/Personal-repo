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

def wavelength_to_rgb(nm) -> Tuple:

    '''
    Takes a wavelength of visible light
    between 380 and 780 nanometres inclusive.
    Values outside this range will raise a ValueError.
    Returns a list of corresponding RGB values.
    Based on Dan Bruton's Fortran implementation.
    '''

    # raise error if nm is outside a (very generous)
    # range of visible light
    if nm < 380 or nm > 780:

        raise ValueError("nm argument must be between 380 and 780")

    # a few variables for later use
    gamma = 0.8
    max_intensity = 255
    factor = 0

    # a dictionary with values initialised to 0
    # which will be calculated later
    rgb = {"R": 0, "G": 0, "B": 0}

    # set RGB values depending on ranges of wavelength
    if 380 <= nm <= 439:
        rgb["R"] = -(nm - 440) / (440 - 380)
        rgb["G"] = 0.0
        rgb["B"] = 1.0
    elif 440 <= nm <= 489:
        rgb["R"] = 0.0
        rgb["G"] = (nm - 440) / (490 - 440)
        rgb["B"] = 1.0
    elif 490 <= nm <= 509:
        rgb["R"] = 0.0
        rgb["G"] = 1.0
        rgb["B"] = -(nm - 510) / (510 - 490)
    elif 510 <= nm <= 579:
        rgb["R"] = (nm - 510) / (580 - 510)
        rgb["G"] = 1.0
        rgb["B"] = 0.0
    elif 580 <= nm <= 644:
        rgb["R"] = 1.0
        rgb["G"] = -(nm - 645) / (645 - 580)
        rgb["B"] = 0.0
    elif 645 <= nm <= 780:
        rgb["R"] = 1.0
        rgb["G"] = 0.0
        rgb["B"] = 0.0

    # calculate a factor (value to multiply by)
    # depending on range of nm
    if 380 <= nm <= 419:
        factor = 0.3 + 0.7 * (nm - 380) / (420 - 380)
    elif 420 <= nm <= 700:
        factor = 1.0
    elif 701 <= nm <= 780:
        factor = 0.3 + 0.7 * (780 - nm) / (780 - 700)

    # adjust RGB values using various variables if > 0
    # else set to 0
    if rgb["R"] > 0:
        rgb["R"] = int(max_intensity * ((rgb["R"] * factor) ** gamma))
    else:
        rgb["R"] = 0

    if rgb["G"] > 0:
        rgb["G"] = int(max_intensity * ((rgb["G"] * factor) ** gamma))
    else:
        rgb["G"] = 0

    if rgb["B"] > 0:
        rgb["B"] = int(max_intensity * ((rgb["B"] * factor) ** gamma))
    else:
        rgb["B"] = 0

    return (rgb["R"], rgb["G"], rgb["B"])
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

landa=640*10**-6
p=500/(10**-3)
D=np.array([39.5,49.5,59.5,69.5,79.5,89.5,99.5])
Xmin1=np.array([0.233/2,0.350/2,0.439/2,0.541/2,0.636/2,0.681/2,0.741/2])
#tabla1=np.array([[39.5,0.117,1.834],
                 #[49.5,0.175,2.460],
                 #[59.5,0.220,2.936],
                 #[69.5,0.271,3.477],
                 #[79.5,0.318,3.833],
                 #[89.5,0.341,4.432],
                 #[99.5,0.371,4.893]])
#print(tabulate.tabulate(tabla1,["Distancia","Posmin","Posmax"],tablefmt="latex"))
Xmax1=np.array([1.834,2.460,2.936,3.477,3.833,4.432,4.893])
#Xmax2=np.array([0.01,0.130,0.391,0.444,0.489,0.523,0.735])
theta=np.array([324,337,337.8,338.6,341,343.8])-324
LE=np.array(["Naranja","Violeta","Azul","Verde","Amarillo","Rojo"])
errorX=np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05])
errorY=np.array([0.001,0.001,0.001,0.001,0.001,0.001,0.001])
f1,m1,b1,r21,m1val=ajuste_linear(D,Xmin1,errorX,errorY,"D","Xmin","Ajuste linear 1")
print('m1=',m1val,'2a=',landa/m1val)
plt.show()
f2,m2,b2,r22,m2val=ajuste_linear(D,Xmax1,errorX,errorY,"D","Xmax","Ajuste linear 2")
print('m2=',m2val,'p=',m2val/landa)
plt.show()
m1=symbols('m1')
m2=symbols('m2')
theta1=symbols('symbols')
print(properrores([m1],landa/m1,[0.00025449905777460806],np.array([[0.004255357142857135]])))
print('sep')
properrores([m2],m2/landa,[0.001],np.array([[0.05]]))
print(properrores([theta1],theta1/p,[np.pi/180],np.pi/180*np.array([[13,13.8,14.6,17,19.8]])))
fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.3)
axbox = fig.add_axes([0.33, 0.01, 0.1, 0.05])
axslider = fig.add_axes([0.7,0.01,0.24,0.05])
slider = Slider(axslider, "Longitud de onda[nm]",380,750)
text_box = TextBox(axbox, "Anchura de la rendija [mm]", textalignment="center")
x=symbols('x')
I0=1
def Intensidad(val):
    lamb = slider.val * 10 ** -9
    a=float(text_box.text)*10**-3
    b=(np.pi*a*sin(x))/lamb
    rgb=wavelength_to_rgb(slider.val)
    col=(1/255 * rgb[0], 1/255* rgb[1], 1/255 * rgb[2])
    plt=plot(I0 * (sin(b) / b) ** 2,(x,-0.03,0.03),fig=fig,ax=ax,line_color=col,xlabel='x [m]' ,ylabel='I(x) [W/m^2]')
    plt.show()
def clear(a):
    ax.clear()
slider.on_changed(clear)
slider.on_changed(Intensidad)
text_box.on_submit(clear)
text_box.on_submit(Intensidad)
text_box.set_val("0.175")
landaexp=np.array([theta,theta*np.pi/180/p*10**9])
print('longitud de onda experimental=',landaexp)

color=[]
for i in range(1,len(landaexp[1])):
    c0=wavelength_to_rgb(round(landaexp[1][i]))
    c=1/255 * c0[0], 1/255* c0[1], 1/255 * c0[2]
    color.append(c)
x2= Symbol('x2')
plot1=plot_implicit(Eq(x2, 0),(x2,-1,23),show=False,title='Espectro experimental',xlabel='theta',line_color=color[4])
plot2=plot_implicit(Eq(x2, 13),(x2,-1,23),show=False,line_color=color[0])
plot3=plot_implicit(Eq(x2, 13.8),(x2,-1,23),show=False,line_color=color[1])
plot4=plot_implicit(Eq(x2, 14.6),(x2,-1,23),show=False,line_color=color[2])
plot5=plot_implicit(Eq(x2, 17),(x2,-1,23),show=False,line_color=color[3])
plot6=plot_implicit(Eq(x2, 19.8),(x2,-1,23),show=False,line_color=color[4])
plot1.append(plot2[0])
plot1.append(plot3[0])
plot1.append(plot4[0])
plot1.append(plot5[0])
plot1.append(plot6[0])
plot1.show()

